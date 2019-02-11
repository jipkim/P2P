type Trade
    tradeID::Int
    kind::Int
    Pes::Float64
    Peb::Float64
    As::Int #Agent - Seller
    Ab::Int #Agent - Buyer
    costNw::Float64
    function Trade(tradeID, kind, Pes, Peb, As, Ab, costNw)
        w = new(tradeID, kind, Pes, Peb, As, Ab, costNw)
        return w
    end
end


type Agent
    agentID::Int
    location::Int
    kind::Int #1-Prosumer 2-Seller 3-Buyer
    pgMin::Float64
    pgMax::Float64
    costFn::Array{Float64,1}
    Pd::Float64
    Qd::Float64
    function Agent(agentID, location, kind, pgMin, pgMax, costFn, Pd, Qd)
        a = new(agentID, location, kind, pgMin, pgMax, costFn, Pd, Qd)
        return a
    end
end
