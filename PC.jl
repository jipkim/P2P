using JuMP, Gurobi, Ipopt
using DataFrames, JLD, CSVFiles
using MatpowerCases

### Include the file with input data functions
include("NetworkDataType.jl")
include("NetworkLoad.jl")
include("P2PTradeType.jl")
include("TradeFunction.jl")
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
gensetU = root
gensetP = setdiff(genset, gensetU)

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
###
trade_scale = 1e-3 # [MWh]

agents = Agent[]
for b in 1:length(buses)
    if isempty(B_gn[b])
        temp_agent = Agent(b, b, 1, 0, 0, [], gam[b]*buses[b].Pd, gam[b]*buses[b].Qd)
        push!(agents,temp_agent)
    else
        for g in B_gn[b]
            if generators[g].gtype != "Root"
                temp_agent = Agent(b, b, 1, generators[g].Pmin, generators[g].Pmax, generators[g].cost, gam[b]*buses[b].Pd, gam[b]*buses[b].Qd)
                push!(agents,temp_agent)
            end
        end
    end
end
# agents = Agent[]
# for b in 1:length(buses)
#     if isempty(B_gn[b])
#         temp_agent = Agent(b, b, 1, 0, 0, [], gam[b]*buses[b].Pd, gam[b]*buses[b].Qd)
#         push!(agents,temp_agent)
#     elseif generators[B_gn[b][1]].gtype != "Root"
#         g = B_gn[b][1]
#         temp_agent = Agent(b, b, 1, generators[g].Pmin, generators[g].Pmax, generators[g].cost, gam[b]*buses[b].Pd, gam[b]*buses[b].Qd)
#         push!(agents,temp_agent)
#     end
# end


# Trade(tradeID, kind, Pes, Peb, As, Ab, costNw)
trades = generateTrade(agents, trade_scale)
trades_old = generateTrade(agents, trade_scale)

trades_mat = Array{Any}(length(agents),length(agents))
trades_mat_buy = Array{Any}(length(agents),length(agents))
trades_mat_sell = Array{Any}(length(agents),length(agents))
for ii in 1:length(agents), jj in 1:length(agents)
    trades_mat[ii,jj] = []
    trades_mat_buy[ii,jj] = []
    trades_mat_sell[ii,jj] = []
end
trades_selected_mat = Array{Any}(length(buses),length(buses))
trades_selected_mat_buy = Array{Any}(length(buses),length(buses))
trades_selected_mat_sell = Array{Any}(length(buses),length(buses))
for ii in 1:length(buses), jj in 1:length(buses)
    trades_selected_mat[ii,jj] = []
    trades_selected_mat_buy[ii,jj] = []
    trades_selected_mat_sell[ii,jj] = []
end

for w in 1:length(trades)
    push!(trades_mat[trades[w].As,trades[w].Ab],trades[w].tradeID)
end

trades_selected = Int[]
trades_rest = Int[]
trades_selected_old = Int[]
trades_rest_old = Int[]
dispatch = zeros(length(buses))
GenCon = zeros(length(buses))
Generation = zeros(length(buses))
Consumption = zeros(length(buses))
GenCon_old = zeros(length(buses))
dispatch_old = zeros(length(buses))
status = Symbol()
DLMP = Float64[]
DLMP_stack = Any[]
dispatch_stack = Any[]
DLMPInfo = DataFrame()
NodeInfo = DataFrame()
LineInfo = DataFrame()
dispatch_peerG = 0.0*zeros(genset)
param_delta = 1.0
tic()
while true
    #for _ in 1:10
    accepted = Array{Int64}[]
    rejected = Array{Int64}[]
    pg = Float64[]
    Revenue_Sell = Float64[]
    Cost_Buy = Float64[]
    Cost_Network = Float64[]
    Cost_DG = Float64[]
    trades_old[:] = trades[:]
    trades_selected_old = trades_selected
    for w in 1:length(trades)
        trades[w].Pes = 0
        trades[w].Peb = 0
    end

    while true
        accepted = Array{Int64}[]
        rejected = Array{Int64}[]
        pg = Float64[]

        gap = 1e-3
        for agentID in 1:length(agents)
            temp_accepted, temp_rejected, temp_pg, Revenue_Sell, Cost_Buy, Cost_Network, Cost_DG = selectTrade(trades, agents, agentID, gap, trade_scale)
            push!(accepted, temp_accepted)
            push!(rejected, temp_rejected)
            push!(pg, temp_pg)
        end

        deltaP = 5*trade_scale #$5/MWh = $5e-3/kWh
        for w in 1:length(trades)
            if (w in accepted[trades[w].Ab]) && ~(w in accepted[trades[w].As])
                if trades[w].Peb > trades[w].Pes
                    trades[w].Pes += deltaP
                else
                    trades[w].Peb += deltaP
                end
            end
        end
        # diff = sum((w in accepted[trades[w].Ab]) && ~(w in accepted[trades[w].As]) for w in 1:length(trades))
        # if diff == 0
        #     break
        # end

        ###diff = sum(abs(trades[w].Pes - trades_old[w].Pes) for w in 1:length(trades)) + sum(abs(trades[w].Peb - trades_old[w].Peb) for w in 1:length(trades))
        if isempty(trades)
            diff = 0
        else
            diff = sum((w in accepted[trades[w].Ab]) && ~(w in accepted[trades[w].As]) for w in 1:length(trades))
        end
        if diff == 0
            break
        end
    end

    trades_selected = Int[]
    trades_rest = Int[]
    for w in 1:length(trades)
        if (w in accepted[trades[w].Ab]) && (w in accepted[trades[w].As])
            push!(trades_selected, w)
        else
            push!(trades_rest, w)
        end
    end


    for i in 1:length(agents)
        dispatch[i] = ceil(agents[i].Pd/trade_scale) # MW to kW conversion
    end

    for w in trades_selected
        dispatch[agents[trades[w].Ab].location] -= 1
        dispatch[agents[trades[w].As].location] += 1
        GenCon[agents[trades[w].Ab].location] -= 1
        GenCon[agents[trades[w].As].location] += 1
        Consumption[agents[trades[w].Ab].location] -= 1
        Generation[agents[trades[w].As].location] += 1
    end

    GenCon = GenCon*trade_scale
    dispatch = dispatch*trade_scale # Kw to MW conversion
    Generation = Generation*trade_scale
    Consumption = Consumption*trade_scale

    #dispatch_peerG = 0.0*zeros(genset)
    for g in gensetP
        dispatch_peerG[g] = Generation[generators[g].location]
    end

    ###### Calculate Network Charge #####
    SMP = generators[root].cost[end-1]
    status, DLMP, pgextra, NodeInfo, LineInfo, DLMPInfo = calculateDLMP(dispatch_peerG, buses, generators, lines, SMP, gensetP, gensetU);


    ###### Update Network Charge #####
    deltaNw = 1*trade_scale # $1/MWh = $1e-3/kWh
    if status == Symbol(:Optimal)
        for w in 1:length(trades)
            if w in trades_selected
                trades[w].costNw = (DLMP[agents[trades[w].Ab].location] - DLMP[agents[trades[w].As].location])/2
            end
        end
    else
        for w in trades_selected
            trades[w].costNw += param_delta*deltaNw/trade_scale
        end
    end
    push!(dispatch_stack, dispatch)
    push!(DLMP_stack, DLMP)
    if GenCon_old == GenCon && status == Symbol(:Optimal)
        break
    else
        dispatch_old[:] = dispatch[:]
        GenCon_old[:] = GenCon[:]
    end
end
simulation_time = toc()


trades_dis = zeros(length(buses),length(buses));
NWcharge_dis = zeros(length(buses),length(buses));
Echarge_dis = zeros(length(buses),length(buses));
temp = Array{Array{Int64,1}}(length(buses),1)
for i in 1:length(buses)
    temp[i] = []
end

for w in trades_selected
    trades_dis[agents[trades[w].Ab].location,agents[trades[w].As].location] -= 1*trade_scale
    trades_dis[agents[trades[w].As].location,agents[trades[w].Ab].location] += 1*trade_scale
    Echarge_dis[agents[trades[w].Ab].location,agents[trades[w].As].location] -= trades[w].Peb
    Echarge_dis[agents[trades[w].As].location,agents[trades[w].Ab].location] += trades[w].Pes
    NWcharge_dis[agents[trades[w].Ab].location,agents[trades[w].As].location] -= round(trades[w].costNw,4)
    NWcharge_dis[agents[trades[w].As].location,agents[trades[w].Ab].location] += round(trades[w].costNw,4)
    push!(temp[agents[trades[w].Ab].location], w)
    push!(temp[agents[trades[w].As].location], w)
end

dispatch_peerG = 0.0*zeros(genset)
for g in gensetP
    dispatch_peerG[g] = Generation[generators[g].location]
end

###### Calculate Network Charge #####
SMP = generators[root].cost[end-1]
status, DLMP, pg, NodeInfo, LineInfo, DLMPInfo = calculateDLMP(dispatch_peerG, buses, generators, lines, SMP, gensetP, gensetU);

#CSV.write("PeerCentric_DLMP.csv", DLMPInfo);
#CSV.write("PeerCentric_NodeInfo.csv", NodeInfo);
#CSV.write("PeerCentric_LineInfo.csv", LineInfo);

price_mat = zeros(length(buses),length(buses))

#############################################
for w in trades_selected
    push!(trades_selected_mat[agents[trades[w].As].location,agents[trades[w].Ab].location],trades[w].tradeID)
end

price_mat = zeros(length(buses),length(buses));
pricemin = zeros(length(buses),length(buses));
pricemax = zeros(length(buses),length(buses));
for ii in 1:length(buses), jj in 1:length(buses)
    temp_mat = []
    for ww in trades_selected_mat[ii,jj]
        push!(temp_mat,trades[ww].Pes)
        pricemin[ii,jj] = minimum(temp_mat)
        pricemin[jj,ii] = -minimum(temp_mat)
        pricemax[ii,jj] = maximum(temp_mat)
        pricemax[jj,ii] = -maximum(temp_mat)
        if pricemin[ii,jj] == pricemax[ii,jj]
            price_mat[ii,jj] = pricemin[ii,jj]
            price_mat[jj,ii] = -pricemin[ii,jj]
        else
            error("pricemin != pricemax")
        end
    end
end

Nprice_mat = zeros(length(buses),length(buses));
for w in trades_selected
    Nprice_mat[agents[trades[w].Ab].location,agents[trades[w].As].location] = 0.5*round(trades[w].costNw,5)
    Nprice_mat[agents[trades[w].As].location,agents[trades[w].Ab].location] = 0.5*round(trades[w].costNw,5)
end


Nprice_mat_min = zeros(length(buses),length(buses));
Nprice_mat_max = zeros(length(buses),length(buses));
for w = 1:length(trades)
    Nprice_mat_min[agents[trades[w].Ab].location,agents[trades[w].As].location] = min(round(trades[w].costNw,5),Nprice_mat_min[agents[trades[w].Ab].location,agents[trades[w].As].location])
    Nprice_mat_min[agents[trades[w].As].location,agents[trades[w].Ab].location] = min(round(trades[w].costNw,5),Nprice_mat_min[agents[trades[w].As].location,agents[trades[w].Ab].location])
    Nprice_mat_max[agents[trades[w].Ab].location,agents[trades[w].As].location] = max(round(trades[w].costNw,5),Nprice_mat_max[agents[trades[w].Ab].location,agents[trades[w].As].location])
    Nprice_mat_max[agents[trades[w].As].location,agents[trades[w].Ab].location] = max(round(trades[w].costNw,5),Nprice_mat_max[agents[trades[w].As].location,agents[trades[w].Ab].location])
end

Eprice_mat = zeros(length(buses),length(buses));
for w in trades_selected
    Eprice_mat[agents[trades[w].Ab].location,agents[trades[w].As].location] = round(trades[w].Pes,5)
    Eprice_mat[agents[trades[w].As].location,agents[trades[w].Ab].location] = round(trades[w].Peb,5)
end


Eprice_mat_min = zeros(length(buses),length(buses));
Eprice_mat_max = zeros(length(buses),length(buses));
for w = 1:length(trades)
    Eprice_mat_min[agents[trades[w].Ab].location,agents[trades[w].As].location] = min(round(trades[w].Pes,5),Eprice_mat_min[agents[trades[w].Ab].location,agents[trades[w].As].location])
    Eprice_mat_min[agents[trades[w].As].location,agents[trades[w].Ab].location] = min(round(trades[w].Peb,5),Eprice_mat_min[agents[trades[w].As].location,agents[trades[w].Ab].location])
    Eprice_mat_max[agents[trades[w].Ab].location,agents[trades[w].As].location] = max(round(trades[w].Pes,5),Eprice_mat_max[agents[trades[w].Ab].location,agents[trades[w].As].location])
    Eprice_mat_max[agents[trades[w].As].location,agents[trades[w].Ab].location] = max(round(trades[w].Peb,5),Eprice_mat_max[agents[trades[w].As].location,agents[trades[w].Ab].location])
end




Peer_revenue = sum(Echarge_dis,2);
Consumer_cost = 0
Generator_revenue = 0
for bb in B_d
    Consumer_cost -= Peer_revenue[bb]
end
for bb in B_g
    Generator_revenue += Peer_revenue[bb]
end

line_marginrate = zeros(length(lines),1)
for ll in 1:length(lines)
    line_marginrate[ll] = round.((lines[ll].u  - sqrt((LineInfo[:fp][ll])^2 + (LineInfo[:fq][ll])^2))/ lines[ll].u,3)
end
line_loading = 1-line_marginrate
bus_voltage = sqrt.(NodeInfo[:v])

v_marginrate = round.(0.1-abs.(1-sqrt.(NodeInfo[:v])),3)/0.1

Loading = DataFrame(line_loading=vec(line_loading))
Vmargin = DataFrame(v_marginrate=vec(v_marginrate))
#####
loading = zeros(length(lines))
for l in 1:length(lines)
   loading[l] = LineInfo[:flow][l] / lines[l].u
end

flow = LineInfo[:flow]
voltage = sqrt.(NodeInfo[:v])

#######
pnm = trades_dis
demand_P2P = zeros(length(buses))
for m in B_d
    demand_P2P[m] = sum(pnm[n,m] for n in B_g)
end

EP_P = zeros(length(busset))
for w in trades_selected
   EP_P[agents[trades[w].Ab].location] += trades[w].Pes
end
EP_U = zeros(length(busset))
for m in B_d
   EP_U[m] = (buses[m].Pd-demand_P2P[m])*DLMP[m]
end

NUC = zeros(length(busset))
for w in trades_selected
   NUC[agents[trades[w].Ab].location] += 0.5*trades[w].costNw
   NUC[agents[trades[w].As].location] += 0.5*trades[w].costNw
end
UR = EP_U + NUC
CP = EP_U + EP_P + NUC

GC = zeros(length(genset))
for g in genset
   GC[g] = generators[g].cost[end-1] * round(pg[g],5)
end
sum_GC_P = sum(GC[g] for g in gensetP)
sum_GC_U = sum(GC[g] for g in gensetU)
sum_DGR = sum(EP_P)
sum_DGP = sum_DGR - sum_GC_P

UP = sum(UR) - sum_GC_U
CashFlow_mat = DataFrame(EP_P=EP_P,EP_U=EP_U,NUC=NUC,UR=UR,CP=CP)
CashFlow_sum = DataFrame(sum_CP=sum(CP), sum_EP_P=sum(EP_P), sum_EP_U=sum(EP_U), sum_NUC=sum(NUC), sum_DGR=sum_DGR, sum_DGP=sum_DGP, UP=UP)

println("partlevel = ""$partlevel")
println("status = ""$status")
println("max loading = ""$(maximum(loading)*100)""%")
println("max voltage = ""$(round(maximum(sqrt.(NodeInfo[:v])),3))")
println("min voltage = ""$(round(minimum(sqrt.(NodeInfo[:v])),3))")
println("max DLMP = ""$(round(maximum(DLMP),3))")
println("min DLMP = ""$(round(minimum(DLMP),3))")
println("simulation time = ""$simulation_time")
method = "PC"
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

#######
resultfile = "$method$(length(buses))""_""$(Int(round(100*partlevel)))""_""$(Int(1000*trade_scale))""kW"".jld"
#@save "$resultfile"



#######################################
trades_selected_n = Array{Array{Int64,1}}(length(buses),1)
trades_n = Array{Array{Int64,1}}(length(buses),1)
for i in 1:length(buses)
    trades_selected_n[i] = []
    trades_n[i] = []
end
for w in trades_selected
    push!(trades_selected_n[trades[w].Ab],w)
    push!(trades_selected_n[trades[w].As],w)
end
for w in 1:length(trades)
    push!(trades_n[trades[w].Ab],w)
    push!(trades_n[trades[w].As],w)
end
Erho = zeros(length(buses),1)
Ecn = zeros(length(buses),1)
Tvol = zeros(length(buses),1)
NUC_total = zeros(length(buses),1)
for i in 1:length(buses)
    if ~isempty(trades_selected_n[i])
        Tvol[i] = length(trades_selected_n[i])
        Ecn[i] = trade_scale*sum(0.5*(DLMP[trades[w].Ab]-DLMP[trades[w].As]) for w in trades_selected_n[i]) / Tvol[i]
        NUC_total[i] = trade_scale*sum(0.5*(DLMP[trades[w].Ab]-DLMP[trades[w].As]) for w in trades_selected_n[i])
        Erho[i] = sum(trades[w].Pes for w in trades_selected_n[i]) / Tvol[i]
    end
end
println("NUC total = ",sum(NUC_total))
println("Tvol = ", sum(Tvol))
#######################################
