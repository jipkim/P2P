function calculateDLMP(dispatch, buses, generators, lines, SMP, gensetP, gensetU)
   lineset = 1:length(lines)
   busset = 1:length(buses)
   genset = 1:length(generators)

   B_g = []
   for g in genset
      push!(B_g, generators[g].location)
   end
   B_d = setdiff(busset,B_g)

   B_gn = Array{Array{Int64}}(length(buses),1)
   for ii in 1:length(buses)
      B_gn[ii] = Int64[]
   end
   for g in genset
      push!(B_gn[generators[g].location], g)
   end
   #########################
   m = Model(solver = IpoptSolver())

   @variable(m, v[b in busset] >= 0) # variable for voltage square, voltage^2
   @variable(m, a[l in lineset] >= 0) # variable for current square, current^2
   @variable(m, fp[l in lineset])
   @variable(m, fq[l in lineset])
   @variable(m, pg[g in genset])
   @variable(m, qg[g in genset])
   @variable(m, OC >= 0)
   @variable(m, OCgen[g in genset] >= 0)
   #@variable(m, pgextra[g in genset]) ##################

   @objective(m, Max, -OC)

   @constraint(m, OperatingCost, OC == sum(OCgen[g] for g in genset))
   #### Polynomial Generation Cost
   #@constraint(m, GenCost[g in genset], OCgen[g] == sum(generators[g].cost[k]*pg[g]^(length(generators[g].cost)-k) for k in 1:length(generators[g].cost)))
   #### Linear Generation Cost
   #@constraint(m, GenCost[g in setdiff(genset,geninter)], OCgen[g] == generators[g].cost[2]*pg[g] + generators[g].cost[3])
   @constraint(m, GenCost[g in genset], OCgen[g] == generators[g].cost[2]*pg[g] + generators[g].cost[3])

   @NLconstraint(m, LineCapFW[l = 1:length(lineset)], (fp[l])^2 + (fq[l])^2 <= ( lines[l].u)^2 )
   #@constraint(m, LineCapFW[l = lineset], norm([fp[l], fq[l]]) <= lines[l].u )

   @NLconstraint(m, LineCapBW[l = 1:length(lineset)], (fp[l] - a[l]*lines[l].r)^2 + (fq[l] - a[l]*lines[l].x)^2 <= ( lines[l].u)^2 )
   #@constraint(m, LineCapBW[l = lineset], norm([fp[l] - a[l]*lines[l].r, (fq[l] - a[l]*lines[l].x)]) <= lines[l].u )

   @constraint(m, BetweenNodes[l = lineset], v[lines[l].fbus] - 2*(lines[l].r*fp[l] + lines[l].x*fq[l]) + a[l]*(lines[l].r^2 + lines[l].x^2) == v[lines[l].tbus])

   @constraint(m, SOCP[l = lineset], ((fp[l])^2 + (fq[l])^2) <= v[lines[l].fbus] * a[l] )
   #@constraint(m, SOCP[l = lineset], norm([2*fp[l], 2*fq[l], v[lines[l].fbus]-a[l]]) <= v[lines[l].fbus] + a[l])

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

   @constraint(m, genPower[g in gensetP], dispatch[g] == pg[g])
   #@constraint(m, genPower[g in gensetP], dispatch[g] + pgextra[g] == pg[g])
   #@constraint(m, extraPower[g in gensetP], pgextra[g]^2 <= (1e-3 * dispatch[g])^2)


   #tic()
   status = solve(m)
   if status != Symbol(:Optimal)
      return status, [], [], [], [], []
   end




   obj_value = getobjectivevalue(m)
   v = getvalue(v)[:]
   a = getvalue(a)[:]
   fp = getvalue(fp)[:]
   fq = getvalue(fq)[:]
   pg = getvalue(pg)[:]
   qg = getvalue(qg)[:]

   OC = getvalue(OC)
   OCgen = getvalue(OCgen)


   pgb = zeros(length(buses))
   qgb = zeros(length(buses))
   for b in 1:length(buses)
      if ~isempty(B_gn[b])
         pgb[b] = sum(pg[g] for g in B_gn[b])
         qgb[b] = sum(qg[g] for g in B_gn[b])
      end
   end

   dual_OC = getdual(OperatingCost)

   dual_LineCapFW = getdual(LineCapFW)
   dual_LineCapBW = getdual(LineCapBW)
   dual_BetweenNodes = getdual(BetweenNodes)
   # dual_SOCP = getdual(SOCP)

   dual_PBalance = getdual(PBalance)
   dual_QBalance = getdual(QBalance)

   dual_Vmax = getdual(Vmax)
   dual_Vmin = getdual(Vmin)

   dual_PGmax = getdual(PGmax)
   dual_PGmin = getdual(PGmin)
   dual_QGmax = getdual(QGmax)
   dual_QGmin = getdual(QGmin)
   A1 = Float64[]
   for l in lineset
       temp = ((fp[l]^2+fq[l]^2)*lines[l].x+a[l]*fq[l]*(lines[l].r^2-lines[l].x^2)-2*a[l]*fp[l]*lines[l].r*lines[l].x)/(
       (fp[l]^2+fq[l]^2)*lines[l].x-a[l]*fq[l]*(lines[l].r^2+lines[l].x^2))
       push!(A1, temp)
   end

   A2 = Float64[]
   for l in lineset
       temp = ((fp[l]^2+fq[l]^2)*lines[l].r-a[l]*fp[l]*(lines[l].r^2+lines[l].x^2))/(
       (fp[l]^2+fq[l]^2)*lines[l].x-a[l]*fq[l]*(lines[l].r^2+lines[l].x^2))
       push!(A2, temp)
   end

   A3 = Float64[]
   for l in lineset
       temp = (-(fp[l]^2+fq[l]^2)*lines[l].r+a[l]*fp[l]*(lines[l].r^2-lines[l].x^2)+2*a[l]*fq[l]*lines[l].r*lines[l].x)/(
       (fp[l]^2+fq[l]^2)*lines[l].x-a[l]*fq[l]*(lines[l].r^2+lines[l].x^2))
       push!(A3, temp)
   end

   A4 = Float64[]
   for l in lineset
       temp = (2*(fq[l]^3*lines[l].r-fp[l]^3*lines[l].x)+2*fp[l]*fq[l]*(fp[l]*lines[l].r-fq[l]*lines[l].x))/(
       (fp[l]^2+fq[l]^2)*lines[l].x-a[l]*fq[l]*(lines[l].r^2+lines[l].x^2))
       push!(A4, temp)
   end

   A5 = Float64[]
   for l in lineset
       temp = (2*(fq[l]^3*lines[l].r-2*fp[l]^3*lines[l].x)+2*fp[l]*fq[l]*(fp[l]*lines[l].r-fq[l]*lines[l].x))/(
       (fp[l]^2+fq[l]^2)*lines[l].x-a[l]*fq[l]*(lines[l].r^2+lines[l].x^2)) +
       (2*a[l]^2*(fq[l]*lines[l].r^3-fp[l]*lines[l].x^3)-4*a[l]*fp[l]*fq[l]*(lines[l].r^2-lines[l].x^2))/(
       (fp[l]^2+fq[l]^2)*lines[l].x-a[l]*fq[l]*(lines[l].r^2+lines[l].x^2)) +
       (4*a[l]*lines[l].r*lines[l].x*(fp[l]^2-fq[l]^2))/(
       (fp[l]^2+fq[l]^2)*lines[l].x-a[l]*fq[l]*(lines[l].r^2+lines[l].x^2)) +
       (-2*a[l]^2*lines[l].r*lines[l].x*(fp[l]*lines[l].r-fq[l]*lines[l].x))/(
       (fp[l]^2+fq[l]^2)*lines[l].x-a[l]*fq[l]*(lines[l].r^2+lines[l].x^2))
       push!(A5, temp)
   end
#=
   DLMP = Array{Float64,1}(length(buses))
   DLMP[1] = SMP
   Eq40 = zeros(length(buses))
   Eq41 = zeros(length(buses))
   Eq42 = zeros(length(buses))
   Eq43 = zeros(length(buses))
   Eq44 = zeros(length(buses))
   for b in 2:length(buses)
       DLMP[b] = A1[b-1]*DLMP[lines[b-1].tbus] + A2[b-1]*dual_QBalance[b] + A3[b-1]*dual_QBalance[lines[b-1].tbus] + A4[b-1]*dual_LineCapFW[b-1] + A5[b-1]*dual_LineCapBW[b-1]
       #DLMP[b] = A1[b-1]*dual_PBalance[lines[b-1].tbus] + A2[b-1]*dual_QBalance[b] + A3[b-1]*dual_QBalance[lines[b-1].tbus] + A4[b-1]*dual_LineCapFW[b-1] + A5[b-1]*dual_LineCapBW[b-1]
       Eq40[b] = A1[b-1]*DLMP[lines[b-1].tbus]
       #Eq40[b] = A1[b-1]*dual_PBalance[lines[b-1].tbus]
       Eq41[b] = A2[b-1]*dual_QBalance[b]
       Eq42[b] = A3[b-1]*dual_QBalance[lines[b-1].tbus]
       Eq43[b] = A4[b-1]*dual_LineCapFW[b-1]
       Eq44[b] = A5[b-1]*dual_LineCapBW[b-1]
   end
=#
   DLMP = Array{Float64,1}(length(buses))
   #DLMP[1] = SMP
   DLMP[1] = dual_PBalance[1]
   Eq40 = zeros(length(buses))
   Eq41 = zeros(length(buses))
   Eq42 = zeros(length(buses))
   Eq43 = zeros(length(buses))
   Eq44 = zeros(length(buses))
   for l in lineset
       DLMP[lines[l].fbus] = A1[l]*dual_PBalance[lines[l].tbus] + A2[l]*dual_QBalance[lines[l].fbus] + A3[l]*dual_QBalance[lines[l].tbus] + A4[l]*dual_LineCapFW[l] + A5[l]*dual_LineCapBW[l]
       Eq40[lines[l].fbus] = A1[l]*dual_PBalance[lines[l].tbus]
       Eq41[lines[l].fbus] = A2[l]*dual_QBalance[lines[l].fbus]
       Eq42[lines[l].fbus] = A3[l]*dual_QBalance[lines[l].tbus]
       Eq43[lines[l].fbus] = A4[l]*dual_LineCapFW[l]
       Eq44[lines[l].fbus] = A5[l]*dual_LineCapBW[l]
   end

   println()
   println("::Network Info::")
   NodeInfo = DataFrame(pg=round.(pgb,3),qg=round.(qgb,3),v=round.(v,3))
   LineInfo = DataFrame(flow=round.(sqrt.(fp.^2+fq.^2),3),i=round.(a,3),fp=round.(fp,3),fq=round.(fq,3))
   display(NodeInfo)
   println()
   println()
   display(LineInfo)


   #println()
   #println("::Network Info::")
   #NetworkInfo = DataFrame(pg=round.(pg,3),qg=round.(qg,3),flow=round.(fp.^2+fq.^2,3),v=round.(v,3),i=round.(a,3),fp=round.(fp,3),fq=round.(fq,3))
   #display(NetworkInfo)
   #println()
   #println()
   println("::DLMP Info::")
   DLMPInfo = DataFrame(λ_i=round.(DLMP,2), Eq40=round.(Eq40,2), Eq41=round.(Eq41,2), Eq42=round.(Eq42,2), Eq43=round.(Eq43,2), Eq44=round.(Eq44,2))
   display(DLMPInfo)
   println()


   return status, DLMP, pg, NodeInfo, LineInfo, DLMPInfo
end
