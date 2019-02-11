using JuMP, Gurobi, Ipopt
using DataFrames, JLD, CSVFiles
using MatpowerCases

### Include the file with input data functions
include("NetworkDataType.jl")
include("NetworkLoad.jl")
include("DLMP.jl")

testsystem = "AP15busDN"
#testsystem = "case141_OPF"

buses, lines, generators, datamat = NetworkLoad(testsystem)
lineset = 1:length(lines)
busset = 1:length(buses)
genset = 1:length(generators)
buses[1].Vmax = 1
buses[1].Vmin = 1

windset = Int[]
for g in genset
   if datamat["gen"][g,2] == "Wind"
      push!(windset, generators[g].gindex)
   end
end


root = -1
for g in genset
   if generators[g].gtype=="Root"
      root = g
   end
end
if root == -1
   display("No generator at root node")
   generators[1].gtype = "Root"
   root = 1
end

temp = Generator
for g in genset
   if generators[g].gtype == "Root"
      temp = Generator(length(genset)+1, "DG", generators[g].location, generators[g].Pg, generators[g].Qg, generators[g].Qmax,
         generators[g].Qmin, generators[g].Vg, generators[g].mBase, generators[g].status, generators[g].Pmax, generators[g].Pmin,
         generators[g].cost, generators[g].SUcost, generators[g].SDcost, generators[g].RU, generators[g].RD, generators[g].UPtime, generators[g].DNtime)
      #push!(generators, temp)
   end
end
genset = 1:length(generators)

########### define Utility generator set & Peer generator set #########
gensetU = root
gensetP = setdiff(genset, gensetU)
#######################################################################

B_g = []
for g in genset
   if ~any(x->x==generators[g].location, B_g)
      push!(B_g, generators[g].location)
      buses[generators[g].location].Pd = 0
   end
end
B_d = setdiff(busset,B_g)

B_gn = Array{Array{Int64}}(length(buses),1)
for ii in 1:length(buses)
   B_gn[ii] = Int64[]
end
for g in genset
   push!(B_gn[generators[g].location], g)
end
############
partlevel = 0.0
gam = partlevel*ones(busset)
############
m = Model(solver = GurobiSolver())

@variable(m, v[b in busset] >= 0) # variable for voltage square, voltage^2
@variable(m, a[l in lineset] >= 0) # variable for current square, current^2
@variable(m, fp[l in lineset])
@variable(m, fq[l in lineset])
@variable(m, pg[g in genset])
@variable(m, qg[g in genset])
@variable(m, OC >= 0)
@variable(m, OCgen[g in genset] >= 0)
# P2P variable
@variable(m, p[n in busset])
@variable(m, pnm[n in busset,k in busset])


@objective(m, Max, -OC)


# P2P constraints
@constraint(m, BiDirection[n in busset, k in busset], pnm[n,k] + pnm[k,n] == 0)
@constraint(m, SignG[n in B_g, k in busset], pnm[n,k] >= 0)
@constraint(m, SignD[n in B_d, k in busset], pnm[n,k] <= 0)
@constraint(m, NodalBalance[n in busset], p[n] == sum(pnm[n,k] for k in busset) )

@constraint(m, genPower[n in B_g], sum(pg[g] for g in setdiff(B_gn[n], gensetU)) == p[n])
@constraint(m, demand[n in B_d], - gam[n] * buses[n].Pd == p[n])

@constraint(m, OperatingCost, OC == sum(OCgen[g] for g in genset))

#### Polynomial Generation Cost
#@constraint(m, GenCost[g in genset], OCgen[g] == sum(generators[g].cost[k]*pg[g]^(length(generators[g].cost)-k) for k in 1:length(generators[g].cost)))
#### Linear Generation Cost
@constraint(m, GenCost[g in genset], OCgen[g] == generators[g].cost[2]*pg[g] + generators[g].cost[3])

#@constraint(m, LineCapFW[l = lineset], (fp[l])^2 + (fq[l])^2 <= ( lines[l].u)^2 )
@constraint(m, LineCapFW[l = 1:length(lineset)], norm([fp[l], fq[l]]) <= lines[l].u )

#@constraint(m, LineCapBW[l = lineset], (fp[l] - a[l]*lines[l].r)^2 + (fq[l] - a[l]*lines[l].x)^2 <= ( lines[l].u)^2 )
@constraint(m, LineCapBW[l = 1:length(lineset)], norm([fp[l] - a[l]*lines[l].r, (fq[l] - a[l]*lines[l].x)]) <= lines[l].u )

@constraint(m, BetweenNodes[l = lineset], v[lines[l].fbus] - 2*(lines[l].r*fp[l] + lines[l].x*fq[l]) + a[l]*(lines[l].r^2 + lines[l].x^2) == v[lines[l].tbus])

#@constraint(m, SOCP[l = lineset], ((fp[l])^2 + (fq[l])^2) <= v[lines[l].fbus] * a[l] )
@constraint(m, SOCP[l = lineset], norm([2*fp[l], 2*fq[l], v[lines[l].fbus]-a[l]]) <= v[lines[l].fbus] + a[l])

if testsystem == "AP15busDN" || testsystem == "IEEE33busDN" || testsystem == "ISONE8busTN" || testsystem == "IEEE33busDN2"
   @constraint(m, v[1] == 1) #Voltage constraint for root node
end

@constraint(m, PBalance[b = busset], sum(fp[l] for l in buses[b].outline)  - sum(fp[l] - lines[l].r*a[l] for l in buses[b].inline) - sum(pg[g] for g in B_gn[b]) + buses[b].Pd + v[b]*buses[b].Gs == 0)
@constraint(m, QBalance[b = busset], sum(fq[l] for l in buses[b].outline)  - sum(fq[l] - lines[l].x*a[l] for l in buses[b].inline) - sum(qg[g] for g in B_gn[b])+ buses[b].Qd - v[b]*buses[b].Bs == 0)

@constraint(m, Vmax[b = busset], v[b] <= buses[b].Vmax^2) #upper limit constraint for voltage square
@constraint(m, Vmin[b = busset], v[b] >= buses[b].Vmin^2) #lower limit constraint for voltage square


@constraint(m, PGmin[g = genset], pg[g] >= generators[g].Pmin)
@constraint(m, PGmax[g = genset], pg[g] <= generators[g].Pmax)
@constraint(m, QGmin[g = genset], qg[g] >= generators[g].Qmin)
@constraint(m, QGmax[g = genset], qg[g] <= generators[g].Qmax)


#print(m)
show(m)
println()
tic()
status = solve(m)
show(m)
println()
simulation_time = toc()
println()
println("Julia v",VERSION)
println("JuMP v",Pkg.installed("JuMP"))
println("Gurobi v",Gurobi.getlibversion())
#print("CPLEX v",CPLEX.version())


obj_value = getobjectivevalue(m)
v = getvalue(v)[:]
a = getvalue(a)[:]
fp = getvalue(fp)[:]
fq = getvalue(fq)[:]
pg = getvalue(pg)[:]
qg = getvalue(qg)[:]



OC = getvalue(OC)
OCgen = getvalue(OCgen)[:]


pgb = zeros(length(buses))
qgb = zeros(length(buses))
for b in 1:length(buses)
   if ~isempty(B_gn[b])
      pgb[b] = sum(pg[g] for g in B_gn[b])
      qgb[b] = sum(qg[g] for g in B_gn[b])
   end
end
pnm = getvalue(pnm)[:,:]
pnm[(abs.(pnm).<1e-3)]=0
result_pnm = pnm


DLMPInfo = DataFrame()
NodeInfo = DataFrame()
LineInfo = DataFrame()

dispatch = pg
SMP = generators[root].cost[end-1]
status1, DLMP, pg, NodeInfo, LineInfo, DLMPInfo = calculateDLMP(dispatch, buses, generators, lines, SMP, gensetP, gensetU);

cn = zeros(length(buses), length(buses))
for n in 1:length(buses), k in 1:length(buses)
    cn[n,k] = (DLMP[k]-DLMP[n])/2
end


line_marginrate = zeros(length(lines),1)
for ll in 1:length(lines)
    line_marginrate[ll] = round.((lines[ll].u  - sqrt((fp[ll])^2 + (fq[ll])^2))/ lines[ll].u,3)
end
line_loading = 1-line_marginrate
bus_voltage = sqrt.(v)

v_marginrate = round.(0.1-abs.(1-sqrt.(v)),3)/0.1

Loading = DataFrame(line_loading=vec(line_loading))
Vmargin = DataFrame(v_marginrate=vec(v_marginrate))

loading = zeros(length(lines))
for l in 1:length(lines)
   loading[l] = LineInfo[:flow][l] / lines[l].u
end

flow = LineInfo[:flow]
voltage = sqrt.(NodeInfo[:v])
##########
EP_P = zeros(length(busset))
for m in B_d
   EP_P[m] = sum(sum(generators[g].cost[end-1]*pnm[n,m] for g in setdiff(B_gn[n],gensetU)) for n in B_g if ~isempty(setdiff(B_gn[n],gensetU)))
end
EP_U = zeros(length(busset))
for m in B_d
   EP_U[m] = (1-gam[m])*buses[m].Pd*DLMP[m]
end
NUC = zeros(length(busset))
for m in B_d
   NUC[m] = sum(pnm[n,m]*(DLMP[m]-DLMP[n]) for n in B_g)
end
UR = EP_U + NUC
CP = EP_U + EP_P + NUC

GC = OCgen
sum_GC_P = sum(GC[g] for g in gensetP)
sum_GC_U = sum(GC[g] for g in gensetU)
sum_DGR = sum(GC[g] for g in gensetP)
sum_DGP = sum_DGR - sum_GC_P
UP = sum(UR) - sum_GC_U

UP = sum(UR) - sum_GC_U
CashFlow_mat = DataFrame(EP_P=EP_P,EP_U=EP_U,NUC=NUC,UR=UR,CP=CP)
CashFlow_sum = DataFrame(sum_CP=sum(CP), sum_EP_P=sum(EP_P), sum_EP_U=sum(EP_U), sum_NUC=sum(NUC), sum_DGR=sum_DGR, sum_DGP=sum_DGP, UP=UP)

println("partlevel = ""$partlevel")
println("status = ""$status")
println("max loading = ""$(maximum(loading)*100)""%")
println("max voltage = ""$(round(maximum(sqrt.(v)),3))")
println("min voltage = ""$(round(minimum(sqrt.(v)),3))")
println("max DLMP = ""$(round(maximum(DLMP),3))")
println("min DLMP = ""$(round(minimum(DLMP),3))")

method = "SC"
voltagefile = "$method$(length(buses))""_""$(Int(100*partlevel))""_voltage"".csv"
loadingfile = "$method$(length(buses))""_""$(Int(100*partlevel))""_loading"".csv"
DLMPfile = "$method$(length(buses))""_""$(Int(100*partlevel))""_DLMP"".csv"
cashflowmatfile = "$method$(length(buses))""_""$(Int(100*partlevel))""_cfmat"".csv"
cashflowsumfile = "$method$(length(buses))""_""$(Int(100*partlevel))""_cfsum"".csv"
# CSVFiles.save(voltagefile, DataFrame(voltage=voltage))
# CSVFiles.save(loadingfile, DataFrame(loading=100*loading))
# CSVFiles.save(DLMPfile, DataFrame(DLMP=DLMP))
# CSVFiles.save(cashflowmatfile, CashFlow_mat)
# CSVFiles.save(cashflowsumfile, CashFlow_sum)

# CSVFiles.save("LineInfo.csv",LineInfo)
##########
#resultfile = "output/$method$(length(buses))""_""$(Int(100*partlevel))"".jld"
resultfile = "$method$(length(buses))""_""$(Int(100*partlevel))"".jld"
# @save "$resultfile"

#######################################
Erho = zeros(length(buses),1)
Ecn = zeros(length(buses),1)

Tvol = sqrt.(sum(pnm,2).^2)
for n in 1:length(busset)
   if Tvol[n] != 0
      Ecn[n] = sum((DLMP[m]-DLMP[n])*pnm[n,m]/2 for m in 1:length(busset)) / Tvol[n]
   end
end
for n in 1:length(busset)
   if Tvol[n] != 0
      Erho[n] = sum((DLMP[m]+DLMP[n])*pnm[n,m]/2 for m in 1:length(busset)) / sum(pnm[n,m] for m in 1:length(busset))
   end
end

NUC_total = zeros(length(buses),1)
for n in 1:length(busset)
   if Tvol[n] != 0
      NUC_total[n] = sum((DLMP[m]-DLMP[n])*pnm[n,m]/2 for m in 1:length(busset))
   end
end
println("NUC total = ",sum(NUC_total))
println("Tvol = ", sum(Tvol))
