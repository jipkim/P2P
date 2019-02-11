function generateTrade(agents::Array{Agent,1}, trade_scale)
   trades = Trade[]
   tradeID = 1
   for s in 1:length(agents), b in 1:length(agents)
      #if s!=b
         for k in 1:min(floor(agents[s].pgMax/trade_scale), ceil(agents[b].Pd/trade_scale))
            w = Trade(tradeID, 1, 0.0, 0.0, s, b, 0)
            push!(trades, w)
            tradeID+=1
         end
      #end
   end
   return trades
end


function selectTrade(W, A, agentID, gap, trade_scale)
    # Trade Set W
    # agaent ∃a ∈ A
    m = Model(solver = GurobiSolver(MIPGap=gap,OutputFlag=0))

    @variable(m, u[1:length(W)], Bin) # whether to take the trade ω or not
    @variable(m, pg >= 0)
    @variable(m, Revenue_Sell)
    @variable(m, Cost_Buy)
    @variable(m, Cost_Network)
    @variable(m, Cost_DG)

    @constraint(m, Revenue_Sell == sum(W[i].Pes * u[i] for i in 1:length(W) if W[i].As == agentID))
    @constraint(m, Cost_Buy == sum(W[i].Peb * u[i] for i in 1:length(W) if W[i].Ab == agentID))
    @constraint(m, Cost_Network == trade_scale*sum(W[i].costNw * u[i] for i in 1:length(W) if W[i].Ab == agentID) + trade_scale*sum(W[i].costNw * u[i] for i in 1:length(W) if W[i].As == agentID))
    #@constraint(m, Cost_Network == 0.5*sum(W[i].costNw * u[i] for i in 1:length(W) if W[i].Ab == agentID) + 0.5*sum(W[i].costNw * u[i] for i in 1:length(W) if W[i].As == agentID))
    #@constraint(m, Cost_Network == sum(W[i].costNw * u[i] for i in 1:length(W) if W[i].Ab == agentID) + sum(W[i].costNw * u[i] for i in 1:length(W) if W[i].As == agentID))
    #@constraint(m, Cost_Network == sum(W[i].costNw * u[i] for i in 1:length(W) if W[i].Ab == agentID) - sum(W[i].costNw * u[i] for i in 1:length(W) if W[i].As == agentID))
    #@constraint(m, Cost_Network == 0.5*sum(W[i].costNw * u[i] for i in 1:length(W) if W[i].As == agentID))
    #@constraint(m, Cost_Network == 0.5*sum(W[i].costNw * u[i] for i in 1:length(W) if W[i].Ab == agentID))
    @constraint(m, Cost_DG == sum(A[agentID].costFn[z] * (pg)^(z-1) for z in 1:length(A[agentID].costFn)))

    @objective(m, Max, Revenue_Sell - Cost_Buy - Cost_DG - Cost_Network)
    #@constraint(m, Generation, pg*1e3 - ceil(A[agentID].Pd*1e3) == sum(u[i] for i in 1:length(W) if W[i].As == agentID)
         #- sum(u[i] for i in 1:length(W) if W[i].Ab == agentID))
    @constraint(m, Generation, pg/trade_scale == sum(u[i] for i in 1:length(W) if W[i].As == agentID))
    #@constraint(m, Consumption, ceil(A[agentID].Pd/trade_scale) == sum(u[i] for i in 1:length(W) if W[i].Ab == agentID))
    @constraint(m, Consumption, floor(A[agentID].Pd/trade_scale) == sum(u[i] for i in 1:length(W) if W[i].Ab == agentID))

    @constraint(m, MinGen, pg >= A[agentID].pgMin)
    @constraint(m, MaxGen, pg <= A[agentID].pgMax)
    for i in 1:length(W)
        if ~(W[i].As == agentID || W[i].Ab == agentID)
            @constraint(m, u[i] == 0)
        end
    end
    #@constraint(m, Only[i in 1:length(W)], u[i] == 0 if ~(W[i].As == agentID || W[i].Ab == agentID) )
    status = solve(m)
    accepted = find(getvalue(u).==1)
    rejected = find(getvalue(u).==0)
    result_pg = getvalue(pg)
    Revenue_Sell = getvalue(Revenue_Sell)
    Cost_Buy = getvalue(Cost_Buy)
    Cost_Network = getvalue(Cost_Network)
    Cost_DG = getvalue(Cost_DG)
    return accepted, rejected, result_pg, Revenue_Sell, Cost_Buy, Cost_Network, Cost_DG
end
