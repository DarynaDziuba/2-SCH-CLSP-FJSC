module TwoSCHMod

# include packages
using DelimitedFiles
using Statistics
using DataFrames       
using CSV
using Random
using PrettyTables
using XLSX
using Profile
using JuMP, Gurobi
using BenchmarkTools

# define tolerance
const MyAbsTol = 1e-12
const MyRelTol = 0.01

# mutable struct for holding CLSP-JS data
# has own constructor CLSPJSData (#periods, #items, #machines)
mutable struct CLSPJSData
    name::String            # name of instance
    T::Int64                # number of periods
    P::Int64                # number of items
    No::Int64               # number of operations
    M::Int64                # number of machines
    D::Matrix{Float64}      # demand matrix D_it
    D_Cap::Matrix{Float64}  # demand matrix in capacity units
    h_orig::Vector{Float64} # original holding cost h_orig_i
    h::Vector{Float64}      # holding cost h_i
    s::Vector{Float64}      # setup cost s_i
    pc::Vector{Float64}     # production cost pc_i
    C::Vector{Vector{Float64}}              # capacity limit C_mt
    C_factors::Matrix{Float64}              # capacity factors
    cumulativeC::Vector{Vector{Float64}}    # cumulative capacity limit C_mt
    I0::Vector{Float64}                     # initial inventory for each item I0_i
    Route::Vector{Vector{Vector{Int64}}}    # item route

    lo::Vector{Int64}       # period of each operation
    io::Vector{Int64}       # item of each operation
    iolo::Matrix{Int64}     # item and period of each operation
    po::Vector{Float64}     # processing time of each operation
    so::Vector{Float64}     # setup time of each operation
    first_operation::Vector{Bool}        # is this the first operation in the item route
    last_operation::Vector{Bool}         # is this the last operation in the item route
    operation_to_machine::Vector{Int64}  # mapping of operation to machine

    # data that must be defined by function 'preprocessData!'
    nonzerodemand::Int64    # number of nonzero demand elements (excluding the first period)
    TBO::Vector{Float64}    # TBO for each item
    averageD_Cap::Vector{Float64}   # average demand in capacity units for each item
    averageDemand::Vector{Float64}  # average demand for each item
    averagePeriodDemand::Float64    # average total demand per period
    averageTBO::Float64     # average TBO of all items
    capTightness::Float64   # average capacity tightness (=averagePeriodDemand/mean per period capacity)
    lumpiness::Float64      # percentage of zero entries in the demand matrix
    stdvDemand::Float64     # standard deviation of the demand

    sequenceOnMachine::Vector{Vector{Int}}      # sequence of operations on each machine
    topological_order::Vector{Int64}            # topological order of operations
    topological_order_reversed::Vector{Int64}   # reversed topological order of operations
    predecessors::Vector{Vector{Int64}}         # there can be at most 2 predecessors: 1 in item route, 1 in machine sequence
    successors::Vector{Vector{Int64}}           # there can be at most 2 successors: 1 in item route, 1 in machine sequence
    has_pred::Vector{Bool}                      # has predecessor
    has_succ::Vector{Bool}                      # has successor
    ind_in_top_order::Vector{Int64}             # index in topological order
    ind_in_top_order_reversed::Vector{Int64}    # index in reversed topological order

    CLSPJSData(name::String, T::Int64, P::Int64, No::Int64, M::Int64) = new(name,T,P,No,M,zeros(Float64, P,T),zeros(Float64, P,T),zeros(Float64, P),zeros(Float64, P),
                                                                        zeros(P),zeros(P),[zeros(T) for _ in 1:M],zeros(M,T),[zeros(T) for _ in 1:M],
                                                                        zeros(P),[[zeros(Int,M) for _ in 1:T] for _ in 1:P],
                                                                        zeros(No),zeros(No),zeros(2,No),zeros(No),zeros(No),zeros(No),zeros(No),zeros(No),
                                                                        0, zeros(P), zeros(P), zeros(P), 0, 0, 0, 0, 0, 
                                                                        [Int[] for _ in 1:M], zeros(No), zeros(No), [zeros(Int64, 2) for _ in 1:No], [zeros(Int64, 2) for _ in 1:No], zeros(No), zeros(No), zeros(No))
end

# mutable struct for holding Job Schedule that represents a conjuncive graph for job shop scheduling
# has own constructor JobSchedule(CLSPJSData)
mutable struct JobSchedule
    data::CLSPJSData                        # corresponding instance data
    duration::Vector{Float64}
    earliestStart::Vector{Float64}
    latestFinish::Vector{Float64}
    slack::Vector{Float64}

    duration_storage1::Vector{Float64}      # storage of backup values in job schedule for try_addLot2Period
    earliestStart_storage1::Vector{Float64}
    latestFinish_storage1::Vector{Float64}
    slack_storage1::Vector{Float64}
    
    duration_storage2::Vector{Float64}      # storage of backup values in job schedule for CASE IIIb
    earliestStart_storage2::Vector{Float64}
    latestFinish_storage2::Vector{Float64}
    slack_storage2::Vector{Float64}
    
    duration_storage3::Vector{Float64}      # storage of backup values in job schedule for try_EXTandNEW
    earliestStart_storage3::Vector{Float64}
    latestFinish_storage3::Vector{Float64}
    slack_storage3::Vector{Float64}

    JobSchedule(data::CLSPJSData) = new(data,
                                    zeros(data.No), zeros(data.No), zeros(data.No), zeros(data.No),
                                    zeros(data.No), zeros(data.No), zeros(data.No), zeros(data.No),
                                    zeros(data.No), zeros(data.No), zeros(data.No), zeros(data.No),
                                    zeros(data.No), zeros(data.No), zeros(data.No), zeros(data.No))
end

# mutable struct for holding the solution for CLSP-JS instance
# has own constructor CLSPJSSolution (CLSPJSData)
mutable struct CLSPJSSolution
    data::CLSPJSData            # corresponding instance data
    itemlist::Vector{Int64}     # list of items - sorted
    periodlist::Vector{Int64}   # list of periods - sorted
    Icost::Float64          # total inventory cost
    Scost::Float64          # total setup cost
    Pcost::Float64          # total production cost
    Ocost::Float64          # total tardiness penalty (former overtime cost)
    Tcost::Float64          # total overall cost
    X::Vector{Vector{Float64}}      # matrix of lotsizes X_it
    Y::Vector{Vector{Bool}}         # matrix of setup variables Y_it
    I::Vector{Vector{Float64}}      # matrix of inventory I_it

    X_storage::Vector{Vector{Float64}}      # storage of backup values in solution for CASE IIIb
    Y_storage::Vector{Vector{Bool}}
    I_storage::Vector{Vector{Float64}}

    total_maxamount::Float64                # use in case of recursive calls in rec_check_feasibilityV2
    continue_simulation::Bool               # use in case of recursive calls in rec_check_feasibilityV2
    init_amount::Float64                    # for rec_check_feasibilityV2 with recursive simulation...

    X_storage3::Vector{Vector{Float64}}     # storage of backup values in solution for try_EXTandNEW
    Y_storage3::Vector{Vector{Bool}}
    I_storage3::Vector{Vector{Float64}}

    O::Vector{Float64}      # overtime in each period
    schedule::JobSchedule   # schedule of the operations

    CLSPJSSolution(data::CLSPJSData) = new(data,collect(Int, 1:data.P),collect(Int, 2:data.T),
                                                0,0,0,0,0,[zeros(Float64,data.T) for _ in 1:data.P],[zeros(Bool,data.T) for _ in 1:data.P],[zeros(Float64,data.T) for _ in 1:data.P],
                                                [zeros(Float64,data.T) for _ in 1:data.P],[zeros(Bool,data.T) for _ in 1:data.P],[zeros(Float64,data.T) for _ in 1:data.P],
                                                #0, false, true, 
                                                0, true, 0.0,
                                                [zeros(Float64,data.T) for _ in 1:data.P],[zeros(Bool,data.T) for _ in 1:data.P],[zeros(Float64,data.T) for _ in 1:data.P],
                                                zeros(data.T),JobSchedule(data))
end

# mutable struct for demand list
mutable struct DemandList
    length::Int64   # number of elements in the demand list (= length of vectors)
    itemIndex::Vector{Int64}    # vector of item indices ordered according to the demand list
    periodIndex::Vector{Int64}  # vector of period indices ordered according to the demand list
    amounts::Vector{Float64}    # vector of demand amounts ordered according to the demand list

    DemandList(length::Int64) = new(length, zeros(Int, length), zeros(Int, length), zeros(Float64, length))
end

# mutable struct for 2-SCH handling and version control
mutable struct _2SCH
    #2-SCH handling variables >>> create this in planDemandList
    data::CLSPJSData                  # corresponding instance data
    cost_amounts::Vector{Float64}   # lot sizes (incremental) for item i reported after checking EXT and NEW options
    postponementDL::DemandList      # postponement demand list
    desactivatePostponement::Bool   # postponement routine must be desactivated in case I when we use inventory
    postpone::Bool                  # demand element should be put in the postponement demand list (=1) or not (=0)
    DL::Bool                        # work with main demand list (=0) or postponement demand list (=1)
    usePostDL::Bool                 # postponement demand list must be checked (=1) or not (=0)
    notAllNow::Bool                 # indicates whether it is possible to add the whole demand element to its period (=0) or not (=1)
    keepInPostDL::Bool              # demand element from postponement demand list is postponed again (=1) or not (=0)
    nextPos::Int64                  # the position after the last element in the postponement demand list
    replan::Bool                    # try to replan DL from scratch with reduced holding cost (=1) or not (=0)
    stop::Bool                      # algorithm is stopped (=1) or not (=0)

    #2-SCH versions control  # define a new immutable struct for this >> global >> create outside
    shiftProduction::Bool           # shift of production routine is enabled (=1) or not (=0)
    postponement::Bool              # postponement routine is enabled (=1) or not (=0)
    threshold::Float64              # threshold value
    cost_reduction::Float64         # percentage by which holding cost is reduced during replanning in case IV
    
    _2SCH(data::CLSPJSData, shiftProduction::Bool=true, postponement::Bool=true, threshold::Float64=0.5, cost_reduction::Float64=0.0) = new(data, zeros(data.T), DemandList(data.nonzerodemand), false, 0, 0, 1, false, false, 1, false, false, shiftProduction, postponement, threshold, cost_reduction)
end

# prints complete job schedule as a table
function printJobSchedule(js::JobSchedule)
    # Prepare table data
    table_data = hcat(
            js.duration, 
            js.earliestStart, 
            js.latestFinish, 
            js.slack
     )    
    # Table header
    header = ([ 
        "Duration", 
        "Earliest Start",  
        "Latest Finish", 
        "Slack"])    
    # Print table
    pretty_table(table_data; header=header, max_num_of_rows=-1, show_row_number=true, display_size=(100,200))
end

# prints partial job schedule as a table with operations in topological order
# can provide start and end operation to print only a part of the schedule
function printPartialJobSchedule(js::JobSchedule, data::CLSPJSData, start_op::Int = 1, end_op::Int = data.No)
    # Find indices of start and end operations in topological order
    start_idx = findfirst(==(start_op), data.topological_order)
    end_idx = findfirst(==(end_op), data.topological_order)
    if start_idx === nothing || end_idx === nothing
        println("Start or end operation not found in topological order.")
        return
    elseif start_idx > end_idx
        println("Start operation appears after end operation in topological order.")
        return
    end
    # Slice the topological order between start and end
    partial_order = data.topological_order[start_idx:end_idx]
    # Prepare table data, explicitly include operation numbers
    op_ids = partial_order
    table_data = hcat(
        op_ids,
        js.duration[partial_order],
        js.earliestStart[partial_order],
        js.latestFinish[partial_order],
        js.slack[partial_order],
        data.successors[partial_order]
    )
    # Table header
    header = ["Operation", "Duration", "Earliest Start", "Latest Finish", "Slack", "Successors"]
    # Print the table
    println("Job schedule from operation $start_op to $end_op in topological order:")
    pretty_table(table_data; header=header, max_num_of_rows=-1, display_size=(100, 200))
end

# Kahn's algorithm for topological sorting
function kahn_topological_sort(data::CLSPJSData)::Vector{Int64}    
    # Step 1: Calculate in-degrees for each job
    in_degree = zeros(Int, data.No)
    for i in 1:data.No
        for successor in data.successors[i]
            if successor > 0
                in_degree[successor] += 1
            end
        end
    end    
    # Step 2: Initialize a queue with nodes of in-degree 0
    queue = []
    for i in 1:data.No
        if in_degree[i] == 0
            push!(queue, i)
        end
    end    
    # Step 3: Process nodes in the queue
    topological_order = []
    while !isempty(queue)
        # Remove a node from the front of the queue
        current = popfirst!(queue)
        push!(topological_order, current)        
        # Reduce in-degree for each successor
        for successor in data.successors[current]
            if successor > 0
                in_degree[successor] -= 1
                # If in-degree becomes 0, add it to the queue
                if in_degree[successor] == 0
                    push!(queue, successor)
                end
            end
        end
    end    
    # Step 4: Check if the graph had a cycle
    if length(topological_order) != data.No
        return "cycle"
    end
    return topological_order
end

# update earliest start, latest finish and slack in the schedule
# start with the first operation in the item route for item i and period t in the forward pass
# start with the last operation in the item route for item i and period t in the backward pass
function update_schedule!(solution::CLSPJSSolution, o_ind::Vector{Int64})::CLSPJSSolution
    ind = solution.data.ind_in_top_order[o_ind[1]]
    # Forward pass: calculate earliest start and finish times
    @inbounds for i in @view solution.data.topological_order[ind:end]
        if (solution.data.last_operation[i]) & (solution.data.lo[i] > 1) 
            solution.schedule.earliestStart[i] = solution.data.cumulativeC[solution.data.operation_to_machine[i]][solution.data.lo[i]-1] #sum(solution.data.C[solution.data.operation_to_machine[i]][1:solution.data.lo[i]-1])
        else
            solution.schedule.earliestStart[i] = 0.0
        end
        if solution.data.has_pred[i]
            # predecessors finish
            for pred in solution.data.predecessors[i]
                if pred > 0
                    predFinish = solution.schedule.earliestStart[pred] + solution.schedule.duration[pred]
                    if predFinish > solution.schedule.earliestStart[i]
                        solution.schedule.earliestStart[i] = predFinish
                    end
                end
            end
        end
    end         
    ind = solution.data.ind_in_top_order_reversed[o_ind[end]]
    # Backward pass: calculate latest start and finish times
    @inbounds for i in @view solution.data.topological_order_reversed[ind:end]  
        solution.schedule.latestFinish[i] = solution.data.cumulativeC[solution.data.operation_to_machine[i]][solution.data.lo[i]]
        if solution.data.has_succ[i]
            # successors start
            for succ in solution.data.successors[i]
                if succ > 0
                    succStart = solution.schedule.latestFinish[succ] - solution.schedule.duration[succ]
                    if succStart < solution.schedule.latestFinish[i]
                        solution.schedule.latestFinish[i] = succStart
                    end
                end
            end
        end
    end
    # Calculate slack
    solution.schedule.slack .= @views solution.schedule.latestFinish .- solution.schedule.earliestStart .- solution.schedule.duration
    return solution
end

# update duration of operations in item route for item i and period t
function change_duration!(solution::CLSPJSSolution, i::Int64, t::Int64)::CLSPJSSolution
    ind = solution.data.Route[i][t]
    @inbounds for o in ind
        solution.schedule.duration[o] = solution.X[i][t]*solution.data.po[o]+solution.data.so[o]*solution.Y[i][t]
    end
    solution = update_schedule!(solution, ind) 
    return solution
end

# master function for checking feasibility of adding 'amount' units of item i in period t
# 'mode' determines which feasibility quantity estimation strategy is used; possible values: "Simple", "Improved_Conservative" (default), "Improved_Exact", "Simulation_Simple", "Simulation_Improved_Conservative", "Simulation_Improved_Exact"
# returns the amount if it can be added to the period completely or maxamount that can be added
function check_feasibility(solution::CLSPJSSolution, i::Int64, t::Int64, amount::Float64, mode::String = "Improved_Conservative")::Float64

    if mode == "Simple"
        maxamount = check_feasibility_simple(solution, i, t, amount)

    elseif mode == "Improved_Conservative"
        maxamount = check_feasibility_improved_conservative(solution, i, t, amount)

    elseif mode == "Improved_Exact"
        maxamount = check_feasibility_improved_exact(solution, i, t, amount)

    elseif mode == "Simulation_Simple"
        maxamount = check_feasibility_simulation_simple(solution, i, t, amount, -1)

    elseif mode == "Simulation_Improved_Conservative"
        maxamount = check_feasibility_simulation_improved_conservative(solution, i, t, amount, -1)

    elseif mode == "Simulation_Improved_Exact"
        maxamount = check_feasibility_simulation_improved_exact(solution, i, t, amount, -1)

    else
        error("Unknown mode for feasibility check: $mode")
    end

    if amount < maxamount
        maxamount = amount
    end
    
    return maxamount
end

# feasibility check with simple formula
# returns the amount if it can be added to the period completely or maxamount that can be added
function check_feasibility_simple(solution::CLSPJSSolution, i::Int64, t::Int64, amount::Float64)::Float64
    maxamount = amount
    ind = solution.data.Route[i][t] 
    sum_proc_time = 0.0
    sum_setup_time = 0.0
    @inbounds for o in ind
        sum_proc_time += solution.data.po[o]
        sum_setup_time += solution.data.so[o]
        if solution.Y[i][t] == true
            max_poss_am = solution.schedule.slack[o] / sum_proc_time
        else
            max_poss_am = (solution.schedule.slack[o] - sum_setup_time) / sum_proc_time
        end
        if max_poss_am < maxamount
            maxamount = max_poss_am
        end
    end

    if maxamount < MyAbsTol
        maxamount = 0.0
    end

    return maxamount
end

# feasibility check with improved formula that computes additional slack conservatively
# returns the amount if it can be added to the period completely or maxamount that can be added
function check_feasibility_improved_conservative(solution::CLSPJSSolution, i::Int64, t::Int64, amount::Float64)::Float64
    maxamount = amount
    ind = solution.data.Route[i][t]
    sum_proc_time = zeros(Float64, length(ind))
    sum_setup_time = zeros(Float64, length(ind))
    d = 0.0
    max_poss_am = 1e+24

    @inbounds for pos in eachindex(ind)
        o = ind[pos]        
        # computing sum of procesing times and sum of setup times
        sum_setup_time[pos] = solution.data.so[o]
        sum_proc_time[pos] = solution.data.po[o]
        if pos > 1
            sum_proc_time[pos] += sum_proc_time[pos-1]
            sum_setup_time[pos] += sum_setup_time[pos-1]
        end
        if solution.Y[i][t] == true
            if !solution.data.first_operation[o]
                for pos_prime in 1:pos-1
                    o_prime = ind[pos_prime]
                    sum_proc_time_pos_prime_pos = sum_proc_time[pos] - sum_proc_time[pos_prime]
                    d = max(0, solution.schedule.slack[o_prime] - solution.schedule.slack[o])
                    q = min(d/sum_proc_time[pos_prime], max(0, solution.schedule.slack[o])/sum_proc_time_pos_prime_pos)
                    q += max(0, (solution.schedule.slack[o] - (d/sum_proc_time[pos_prime])*sum_proc_time_pos_prime_pos)/sum_proc_time[pos])
                    if q < max_poss_am
                        max_poss_am = q
                    end
                end
            else
                max_poss_am = solution.schedule.slack[o] / solution.data.po[o]
            end
        else
            if !solution.data.first_operation[o]
                for pos_prime in 1:pos-1
                    o_prime = ind[pos_prime]
                    sum_proc_time_pos_prime_pos = sum_proc_time[pos] - sum_proc_time[pos_prime]
                    sum_setup_time_pos_prime_pos = sum_setup_time[pos] - sum_setup_time[pos_prime]
                    d = max(0, solution.schedule.slack[o_prime] - solution.schedule.slack[o])
                    # Case 1
                    if d <= sum_setup_time[pos_prime]
                        q = max(0, (solution.schedule.slack[o] - sum_setup_time[pos] + d)/sum_proc_time[pos])
                    # Case 2
                    else
                        q1 = (d-sum_setup_time[pos_prime])/sum_proc_time[pos_prime]
                        q2 = max(0,(solution.schedule.slack[o]-sum_setup_time_pos_prime_pos)/sum_proc_time_pos_prime_pos)
                        q = min(q1, q2)
                        q += max(0, (solution.schedule.slack[o] - sum_setup_time_pos_prime_pos - (d-sum_setup_time[pos_prime])*(sum_proc_time_pos_prime_pos/sum_proc_time[pos_prime])))/sum_proc_time[pos]
                    end

                    if q < max_poss_am
                        max_poss_am = q
                    end
                end
            else
                max_poss_am = (solution.schedule.slack[o] - solution.data.so[o])/ solution.data.po[o]
            end
        end      
        if max_poss_am < maxamount
            maxamount = max_poss_am
        end
    end

    if maxamount < MyAbsTol
        maxamount = 0.0
    end

    return maxamount
end

# feasibility check with simulation that uses the simple formula, 
# returns the amount if it can be added to the period completely or maxamount that can be added
function check_feasibility_simulation_simple(solution::CLSPJSSolution, i::Int64, t::Int64, amount::Float64, counter_value::Int64 = -1)::Float64
    maxamount = amount
    ind = solution.data.Route[i][t]      # get the indices of operations (all levels) for item i in period t
    remaining_amount = amount
    counter = 1
    total_maxamount = 0.0
    sim = false
    if counter_value == -1
        counter_value = 1000000
    end

    backupX = solution.X[i][t]
    backupY = solution.Y[i][t]
   
    while (remaining_amount > MyAbsTol)         
        if counter > 1
            # Simulate the change
            solution.X[i][t] += maxamount
            solution.Y[i][t] = true
            solution = change_duration!(solution, i, t)
            sim = true
            # update maxamount with remaining amount after simulation
            maxamount = remaining_amount
        end
        
        sum_proc_time = 0.0
        sum_setup_time = 0.0
   
        @inbounds for o in ind
            sum_proc_time += solution.data.po[o]
            sum_setup_time += solution.data.so[o]
            if solution.Y[i][t] == true
                max_poss_am = solution.schedule.slack[o] / sum_proc_time
            else
                max_poss_am = (solution.schedule.slack[o] - sum_setup_time) / sum_proc_time
            end
            if max_poss_am < maxamount
                maxamount = max_poss_am
            end
        end

        counter += 1

        if maxamount >= MyAbsTol
            remaining_amount -= maxamount
            total_maxamount += maxamount
        end
        
        if maxamount < MyAbsTol
            maxamount = 0.0
            break
        end

        if counter == counter_value
            break
        end

    end

    if sim
        # undo simulation
        solution.X[i][t] = backupX
        solution.Y[i][t] = backupY
        solution = change_duration!(solution, i, t)
    end

    if amount < total_maxamount
        return amount
    else
        return total_maxamount
    end
end

# feasibility check with simulation that uses the improved formula with conservative additional slack
# returns the amount if it can be added to the period completely or maxamount that can be added
function check_feasibility_simulation_improved_conservative(solution::CLSPJSSolution, i::Int64, t::Int64, amount::Float64, counter_value::Int64 = -1)::Float64
    maxamount = amount
    ind = solution.data.Route[i][t]      # get the indices of operations (all levels) for item i in period t
    remaining_amount = amount
    counter = 1
    total_maxamount = 0.0
    sim = false
    if counter_value == -1
        counter_value = 1000000
    end

    backupX = solution.X[i][t]
    backupY = solution.Y[i][t]
   
    while (remaining_amount > MyAbsTol)         
        if counter > 1
            # Simulate the change
            solution.X[i][t] += maxamount
            solution.Y[i][t] = true
            solution = change_duration!(solution, i, t)
            sim = true
            # update maxamount with remaining amount after simulation
            maxamount = remaining_amount
        end
        
        sum_proc_time = zeros(Float64, length(ind))
        sum_setup_time = zeros(Float64, length(ind))
        d = 0.0
        max_poss_am = 1e+24
   
        @inbounds for pos in eachindex(ind)
            o = ind[pos]
            
            # computing sum of procesing times and sum of setup times
            sum_setup_time[pos] = solution.data.so[o]
            sum_proc_time[pos] = solution.data.po[o]
            if pos > 1
                sum_proc_time[pos] += sum_proc_time[pos-1]
                sum_setup_time[pos] += sum_setup_time[pos-1]
            end

            if solution.Y[i][t] == true
                if !solution.data.first_operation[o]
                    for pos_prime in 1:pos-1
                        o_prime = ind[pos_prime]
                        sum_proc_time_pos_prime_pos = sum_proc_time[pos] - sum_proc_time[pos_prime]
                        d = max(0, solution.schedule.slack[o_prime] - solution.schedule.slack[o])
                        q = min(d/sum_proc_time[pos_prime], max(0, solution.schedule.slack[o])/sum_proc_time_pos_prime_pos)
                        q += max(0, (solution.schedule.slack[o] - (d/sum_proc_time[pos_prime])*sum_proc_time_pos_prime_pos)/sum_proc_time[pos])
                        if q < max_poss_am
                            max_poss_am = q
                        end
                    end
                else
                    max_poss_am = solution.schedule.slack[o] / solution.data.po[o]
                end
            else
                if !solution.data.first_operation[o]
                    for pos_prime in 1:pos-1
                        o_prime = ind[pos_prime]
                        sum_proc_time_pos_prime_pos = sum_proc_time[pos] - sum_proc_time[pos_prime]
                        sum_setup_time_pos_prime_pos = sum_setup_time[pos] - sum_setup_time[pos_prime]
                        d = max(0, solution.schedule.slack[o_prime] - solution.schedule.slack[o])
                        # Case 1
                        if d <= sum_setup_time[pos_prime]
                            q = max(0, (solution.schedule.slack[o] - sum_setup_time[pos] + d)/sum_proc_time[pos])
                        # Case 2
                        else
                            q1 = (d-sum_setup_time[pos_prime])/sum_proc_time[pos_prime]
                            q2 = max(0,(solution.schedule.slack[o]-sum_setup_time_pos_prime_pos)/sum_proc_time_pos_prime_pos)
                            q = min(q1, q2)
                            q += max(0, (solution.schedule.slack[o] - sum_setup_time_pos_prime_pos - (d-sum_setup_time[pos_prime])*(sum_proc_time_pos_prime_pos/sum_proc_time[pos_prime])))/sum_proc_time[pos]
                        end
    
                        if q < max_poss_am
                            max_poss_am = q
                        end
    
                    end
                else
                    max_poss_am = (solution.schedule.slack[o] - solution.data.so[o])/ solution.data.po[o]
                end
            end    
          
            if max_poss_am < maxamount
                maxamount = max_poss_am
            end
        end

        counter += 1

        if maxamount >= MyAbsTol
            remaining_amount -= maxamount
            total_maxamount += maxamount
        end
        
        if maxamount < MyAbsTol
            maxamount = 0.0
            break
        end

        if counter == counter_value
            break
        end

    end

    if sim
        # undo simulation
        solution.X[i][t] = backupX
        solution.Y[i][t] = backupY
        solution = change_duration!(solution, i, t)
    end

    if amount < total_maxamount
        return amount
    else
        return total_maxamount
    end
end

# compute matrix d_{i,i} with exact additional slack based on the longest path
function compute_correct_dii(solution::CLSPJSSolution, i::Int64, t::Int64)::Matrix{Float64}
    ind = solution.data.Route[i][t]
    pos1 = solution.data.ind_in_top_order[ind[1]]
    pos_end = solution.data.ind_in_top_order[ind[end]]
    
    tau = 1e+100*ones(Float64,(length(ind), pos_end-pos1+1))
    index_prime = 0

    @inbounds for o in @view solution.data.topological_order[pos1:pos_end-1]
        pos_o = solution.data.ind_in_top_order[o]
        if o == ind[index_prime+1]
            index_prime += 1
            tau[index_prime, pos_o-pos1+1] = 0.0 # set to zero for itself
        end

        if (solution.data.has_succ[o])
            for succ in solution.data.successors[o]
                if (succ > 0)
                    pos_succ = solution.data.ind_in_top_order[succ]
                    if pos_succ <= pos_end
                        for j in 1:index_prime
                            if ((tau[j,pos_succ-pos1+1] < tau[j, pos_o - pos1+1] + solution.schedule.duration[o]) || (tau[j, pos_succ-pos1+1] ==  1e+100)) &&(tau[j, pos_o - pos1+1] < 1e+100)
                                tau[j, pos_succ-pos1+1] = tau[j, pos_o - pos1+1] + solution.schedule.duration[o]
                            end
                        end
                    end
                end
            end
        end
    end

    dii = zeros(Float64, (length(ind), length(ind)))
    for index_prime in 1:length(ind)-1
        for index in index_prime+1:length(ind)
            o_index = solution.data.ind_in_top_order[ind[index]]
            if tau[index_prime, o_index - pos1+1] < 1e+100
                dii[index_prime, index] = solution.schedule.earliestStart[ind[index]] - solution.schedule.earliestStart[ind[index_prime]] - tau[index_prime, o_index-pos1+1]
            end
        end
    end

    return dii
end

# feasibility check with improved formula that computes additional slack exactly
# returns the amount if it can be added to the period completely or maxamount that can be added
function check_feasibility_improved_exact(solution::CLSPJSSolution, i::Int64, t::Int64, amount::Float64)::Float64
    maxamount = amount
    ind = solution.data.Route[i][t]
    sum_proc_time = zeros(Float64, length(ind))
    sum_setup_time = zeros(Float64, length(ind))
    dii = compute_correct_dii(solution, i, t)
    max_poss_am = 1e+24

    @inbounds for pos in eachindex(ind)
        o = ind[pos]
        
        # computing sum of procesing times and sum of setup times
        sum_setup_time[pos] = solution.data.so[o]
        sum_proc_time[pos] = solution.data.po[o]
        if pos > 1
            sum_proc_time[pos] += sum_proc_time[pos-1]
            sum_setup_time[pos] += sum_setup_time[pos-1]
        end

        if solution.Y[i][t] == true
            if !solution.data.first_operation[o]
                for pos_prime in 1:pos-1
                    d = dii[pos_prime, pos]
                    sum_proc_time_pos_prime_pos = sum_proc_time[pos] - sum_proc_time[pos_prime]
                    q = min(d/sum_proc_time[pos_prime], max(0, solution.schedule.slack[o])/sum_proc_time_pos_prime_pos)
                    q += max(0, (solution.schedule.slack[o] - (d/sum_proc_time[pos_prime])*sum_proc_time_pos_prime_pos)/sum_proc_time[pos])
                    if q < max_poss_am
                        max_poss_am = q
                    end
                end
            else
                max_poss_am = solution.schedule.slack[o] / solution.data.po[o]
            end
        else
            if !solution.data.first_operation[o]
                for pos_prime in 1:pos-1
                    d = dii[pos_prime, pos]
                    sum_proc_time_pos_prime_pos = sum_proc_time[pos] - sum_proc_time[pos_prime]
                    sum_setup_time_pos_prime_pos = sum_setup_time[pos] - sum_setup_time[pos_prime]
                    # Case 1
                    if d <= sum_setup_time[pos_prime]
                        q = max(0, (solution.schedule.slack[o] - sum_setup_time[pos] + d)/sum_proc_time[pos])
                    # Case 2
                    else
                        q1 = (d-sum_setup_time[pos_prime])/sum_proc_time[pos_prime]
                        q2 = max(0,(solution.schedule.slack[o]-sum_setup_time_pos_prime_pos)/sum_proc_time_pos_prime_pos)
                        q = min(q1, q2)
                        q += max(0, (solution.schedule.slack[o] - sum_setup_time_pos_prime_pos - (d-sum_setup_time[pos_prime])*(sum_proc_time_pos_prime_pos/sum_proc_time[pos_prime])))/sum_proc_time[pos]
                    end

                    if q < max_poss_am
                        max_poss_am = q
                    end

                end
            else
                max_poss_am = (solution.schedule.slack[o] - solution.data.so[o])/ solution.data.po[o]
            end
        end
      
        if max_poss_am < maxamount
            maxamount = max_poss_am
        end
    end

    if maxamount < MyAbsTol
        maxamount = 0.0
    end

    return maxamount
end

# feasibility check with simulation that uses the improved formula with exact additional slack
# returns the amount if it can be added to the period completely or maxamount that can be added
function check_feasibility_simulation_improved_exact(solution::CLSPJSSolution, i::Int64, t::Int64, amount::Float64, counter_value::Int64 = -1)::Float64
    maxamount = amount
    ind = solution.data.Route[i][t]      # get the indices of operations (all levels) for item i in period t
    remaining_amount = amount
    counter = 1
    total_maxamount = 0.0
    sim = false
    if counter_value == -1
        counter_value = 1000000
    end

    backupX = solution.X[i][t]
    backupY = solution.Y[i][t]
   
    while (remaining_amount > MyAbsTol)        
        if counter > 1
            # Simulate the change
            solution.X[i][t] += maxamount
            solution.Y[i][t] = true
            solution = change_duration!(solution, i, t)
            sim = true
            # update maxamount with remaining amount after simulation
            maxamount = remaining_amount
        end
        
        sum_proc_time = zeros(Float64, length(ind))
        sum_setup_time = zeros(Float64, length(ind))
        d = 0.0
        dii = compute_correct_dii(solution, i, t)
        max_poss_am = 1e+24
   
        @inbounds for pos in eachindex(ind)
            o = ind[pos]
            
            # computing sum of procesing times and sum of setup times
            sum_setup_time[pos] = solution.data.so[o]
            sum_proc_time[pos] = solution.data.po[o]
            if pos > 1
                sum_proc_time[pos] += sum_proc_time[pos-1]
                sum_setup_time[pos] += sum_setup_time[pos-1]
            end

            if solution.Y[i][t] == true
                if !solution.data.first_operation[o]
                    for pos_prime in 1:pos-1
                        sum_proc_time_pos_prime_pos = sum_proc_time[pos] - sum_proc_time[pos_prime]
                        d = dii[pos_prime, pos]
                        q = min(d/sum_proc_time[pos_prime], max(0, solution.schedule.slack[o])/sum_proc_time_pos_prime_pos)
                        q += max(0, (solution.schedule.slack[o] - (d/sum_proc_time[pos_prime])*sum_proc_time_pos_prime_pos)/sum_proc_time[pos])
                        if q < max_poss_am
                            max_poss_am = q
                        end
                    end
                else
                    max_poss_am = solution.schedule.slack[o] / solution.data.po[o]
                end
            else
                if !solution.data.first_operation[o]
                    for pos_prime in 1:pos-1
                        sum_proc_time_pos_prime_pos = sum_proc_time[pos] - sum_proc_time[pos_prime]
                        sum_setup_time_pos_prime_pos = sum_setup_time[pos] - sum_setup_time[pos_prime]
                        d = dii[pos_prime, pos]
                        # Case 1
                        if d <= sum_setup_time[pos_prime]
                            q = max(0, (solution.schedule.slack[o] - sum_setup_time[pos] + d)/sum_proc_time[pos])
                        # Case 2
                        else
                            q1 = (d-sum_setup_time[pos_prime])/sum_proc_time[pos_prime]
                            q2 = max(0,(solution.schedule.slack[o]-sum_setup_time_pos_prime_pos)/sum_proc_time_pos_prime_pos)
                            q = min(q1, q2)
                            q += max(0, (solution.schedule.slack[o] - sum_setup_time_pos_prime_pos - (d-sum_setup_time[pos_prime])*(sum_proc_time_pos_prime_pos/sum_proc_time[pos_prime])))/sum_proc_time[pos]
                        end
    
                        if q < max_poss_am
                            max_poss_am = q
                        end
    
                    end
                else
                    max_poss_am = (solution.schedule.slack[o] - solution.data.so[o])/ solution.data.po[o]
                end
            end    
          
            if max_poss_am < maxamount
                maxamount = max_poss_am
            end
        end

        counter += 1

        if maxamount >= MyAbsTol
            remaining_amount -= maxamount
            total_maxamount += maxamount
        end
        
        if maxamount < MyAbsTol
            maxamount = 0.0
            break
        end

        if counter == counter_value
            break
        end

    end

    if sim
        # undo simulation
        solution.X[i][t] = backupX
        solution.Y[i][t] = backupY
        solution = change_duration!(solution, i, t)
    end

    if amount < total_maxamount
        return amount
    else
        return total_maxamount
    end
end

# plans the first period and updates the job schedule accordingly
function planFirstPeriod!(solution::CLSPJSSolution)::CLSPJSSolution
    # plan demand in the first period and determine the cost of production plan
    solution.Icost = 0.0
    solution.Scost = 0.0
    solution.Pcost = 0.0
    @inbounds for i in 1:solution.data.P
        if (solution.data.D[i,1] > MyAbsTol)
            solution.X[i][1] = solution.data.D[i,1]
            ind = solution.data.Route[i][1]               # get the indices of operations (all levels) for item i in period 1
            @inbounds for o in ind
                solution.schedule.duration[o] = solution.X[i][1]*solution.data.po[o] + solution.data.so[o]
            end
            # update schedule
            solution = update_schedule!(solution, ind)
            @inbounds for o in ind
                m = solution.data.operation_to_machine[o]
                # update tardiness penalty
                solution.Ocost += 10000 * max(0.0, solution.schedule.latestFinish[o] - solution.data.cumulativeC[m][solution.data.lo[o]])
            end
            solution.Y[i][1] = true
            solution.Scost += solution.data.s[i]
            solution.Pcost += solution.data.pc[i] * solution.X[i][1]
        end
    end
    solution.Tcost = solution.Scost + solution.Pcost + solution.Ocost
    solution.Tcost = update_Tcostv2(solution)
    return solution
end

# pre-process data: count nonzero demand elements, compute TBO, etc., setup job schedule, convert capacity to %
function preprocessData!(data::CLSPJSData)::CLSPJSData
    # compute cumulative available capacity
    for m in 1:data.M
        for t in 1:data.T
            data.cumulativeC[m][t] = sum(data.C[m][1:t])
        end
    end
    # count nonzero demand elements
    for i in 1:data.P, t in 2:data.T
        if data.D[i,t] > MyAbsTol
            data.nonzerodemand += 1
        end
    end
    # compute capacity utilization for each demand element
    @inbounds @views for i in 1:data.P, t in 1:data.T
        data.D_Cap[i,t] = sum(data.po[o] for o in data.Route[i][t])*data.D[i,t] / data.C[1][t]
    end
    # compute average demand of each item
    for i in 1:data.P
        for t in 1:data.T
            data.averageDemand[i] += data.D[i,t]
        end
        data.averageDemand[i] /= data.T
    end
    # compute average capacity utilization of each item, i.e., convert demand units to capacity units
    for i in 1:data.P
        for t in 1:data.T
            data.averageD_Cap[i] += data.D_Cap[i,t]
        end
        data.averageD_Cap[i] /= data.T
    end
    # compute average period demand (in capacity units)
    for i in 1:data.P
        for t in 1:data.T
            data.averagePeriodDemand += data.D_Cap[i,t]
        end
        data.averagePeriodDemand /= data.T
    end
    # compute TBO of each item and average TBO
    averageSetupCap = 0
    for i in 1:data.P
        data.TBO[i] = sqrt((2*data.s[i])/(data.h[i]*data.averageDemand[i]))
        averageSetupCap += (sum(data.so[o] for t in 1:data.T for o in data.Route[i][t]) / data.TBO[i])*data.M # multiplied by M to capture the number of operations in item route
    end
    data.averageTBO = mean(data.TBO)
    averageCap = mean(map(sum, data.C)) 
    data.capTightness = (sum(data.averageD_Cap) + averageSetupCap) / averageCap # compute capacity tightness accounting for setup times

    for i in 1:data.P,t in 1:data.T
        if data.D[i,t] <= MyAbsTol
            data.lumpiness += 1
        end
    end
    data.lumpiness /= (data.P*data.T)     # compute lumpiness
    data.stdvDemand = std(data.D)         # compute standard deviation of demand
    return data
end

# read instance from Wolosewich et al. 2015
function W_read_CLSPJSData(instance_path::String)::CLSPJSData
    filename = instance_path*"DonneesGenerales.txt"
    lines = readlines(filename)
    ind = 0
    P = parse(Int64,strip(lines[ind += 1]))
    T = parse(Int64,strip(lines[ind += 1]))
    No = parse(Int64,strip(lines[ind += 1]))
    M = parse(Int64,strip(lines[ind += 1]))
    data = CLSPJSData(instance_path,T,P,No,M)               # create instance data

    filename = instance_path*"ArcsPrecedences.txt"          # read product routes
    lines = readlines(filename) 
    ind = 0
    for t in 1:data.T
        for i in 1:data.P
            data.Route[i][t] .= parse.(Int64, split(lines[ind+=1], " "))
        end
    end

    filename = instance_path*"ArcsSequence.txt"             # read machine sequences
    lines = readlines(filename)
    for m in 1:data.M
        data.sequenceOnMachine[m] = parse.(Int64, split(lines[m], " "))
    end

    filename = instance_path*"DonneesCapacite.txt"          # read capacity utilization parameters
    lines = readlines(filename)
    for m in 1:data.M
        data.C[m] = parse.(Float64, lines)
    end

    filename = instance_path*"DonneesOperation.txt"         # read operation data
    lines = readlines(filename)
    ind = -1
    for o in 1:data.No
        op, data.io[o], data.lo[o], data.operation_to_machine[o], data.po[o], data.so[o] = parse.(Int64, split(lines[ind+=2], " "))
    end

    filename = instance_path*"DonneesProduit.txt"             # read product data
    lines = readlines(filename)
    if length(lines) != data.P
        ind = 0
        for i in 1:data.P
            item, data.pc[i], bc, data.h[i], data.s[i], lead_time = parse.(Float64, split(strip(lines[ind+=1]), " "))
            data.D[i, :] .= parse.(Float64, split(strip(lines[ind+=1]), " "))
        end
    else
        for i in 1:data.P
            item, data.pc[i], bc, data.h[i], data.s[i], lead_time = parse.(Float64, split(lines[i], " ")[1:6])
            data.D[i, :] .= parse.(Float64, split(lines[i], " ")[7:end])
        end
    end
    data.h_orig .= data.h
    data.iolo[1, :] .= data.io
    data.iolo[2, :] .= data.lo 

    # Compute capacity from processing and setup times and demand, apply cap.utilization parameter    
    reqC = zeros(data.T)
    for t in 1:data.T, i in 1:data.P
        for o in data.Route[i][t]
            reqC[t] += data.po[o]*data.D[i,t] + data.so[o]
        end
    end
    for m in 1:data.M
        data.C[m] .*= reqC
    end

    # Determine the first and the last operations
    for i in 1:P, t in 1:T
        data.first_operation[data.Route[i][t][1]] = true
        data.last_operation[data.Route[i][t][end]] = true
    end

    # Define the precedence relations
    for o in 1:data.No
        pred = 1
        # predecessor machine-wise
        pos = findfirst(x -> x == o, data.sequenceOnMachine[data.operation_to_machine[o]])
        if pos !== nothing && pos > 1
            data.predecessors[o][pred] = data.sequenceOnMachine[data.operation_to_machine[o]][pos-1]
            pred += 1
        else
            pred += 1
        end
        # predecessor item-wise
        if !data.first_operation[o]
            pos = findfirst(x -> x == o, data.Route[data.io[o]][data.lo[o]])
            data.predecessors[o][pred] = data.Route[data.io[o],][data.lo[o]][pos-1]
            pred += 1
        end
    end
 
    for o in 1:data.No
        # successors
        if sum(data.predecessors[o]) > 0
            # println("This operation has predecessors")
            data.has_pred[o] = true
            for pred in data.predecessors[o]
                if pred > 0
                    if data.successors[pred][1] == 0
                        data.successors[pred][1] = o
                        data.has_succ[pred] = true
                    elseif data.successors[pred][2] == 0
                        data.successors[pred][2] = o
                        data.has_succ[pred] = true
                    else
                        println("Error with indexing")
                    end
                end
            end
        end
    end

    # Compute the topological order
    data.topological_order = kahn_topological_sort(data)
    data.topological_order_reversed = reverse(data.topological_order)

    data.ind_in_top_order = [findfirst(==(x), data.topological_order) for x in 1:data.No]
    data.ind_in_top_order_reversed = [findfirst(==(x), data.topological_order_reversed) for x in 1:data.No]

    data = preprocessData!(data)
    return data
end

# reset solution: empty production plan and empty schedule, plan the first period
function resetSolution!(solution::CLSPJSSolution)::CLSPJSSolution
    # reset production plan
    @inbounds for p in 1:solution.data.P
        solution.X[p] .= 0.0
        solution.Y[p] .= false
        solution.I[p] .= 0.0
    end
    solution.O .= 0.0

    # reset (set up) operations schedule
    @inbounds for o in 1:solution.data.No
        solution.schedule.duration[o] = 0.0
        m = solution.data.operation_to_machine[o]
        t = solution.data.lo[o]
        if (solution.data.last_operation[o]) & (t>1)    # if this is the last operation in item route
            solution.schedule.earliestStart[o] = solution.data.cumulativeC[m][t-1]
        else
            solution.schedule.earliestStart[o] = 0.0
        end
        solution.schedule.latestFinish[o] =  solution.data.cumulativeC[m][t]
    end
    # Backward pass: calculate latest finish times
    @inbounds for i in @view solution.data.topological_order_reversed[1:end]
        if solution.data.has_succ[i]
            # successors start
            for succ in solution.data.successors[i]
                if succ > 0
                    succStart = solution.schedule.latestFinish[succ] - solution.schedule.duration[succ]
                    if succStart < solution.schedule.latestFinish[i]
                        solution.schedule.latestFinish[i] = succStart
                    end
                end
            end
        end
    end
    solution.schedule.slack .= @views solution.schedule.latestFinish .- solution.schedule.earliestStart .- solution.schedule.duration
    # plan first period
    solution = planFirstPeriod!(solution)
    return solution
end

# create a copy of existing solution
function storeSolution(solution::CLSPJSSolution)::CLSPJSSolution
    solution_copy = CLSPJSSolution(solution.data)
    solution_copy.X = copy(solution.X)
    solution_copy.I = copy(solution.I)
    solution_copy.Y = copy(solution.Y)
    solution_copy.O[:] = copy(solution.O[:])
    solution_copy.Tcost = update_Tcostv2(solution_copy)
    return solution_copy
end

# check if solution satisfies all setup and inventory balance constraints 
function checkSolution(solution::CLSPJSSolution)
    rv = 0
    errmsg = []

    # Convert Vector{Vector{Float64}} to Matrix{Float64} for easier operations
    X_matrix = stack(solution.X, dims=1)
    I_matrix = stack(solution.I, dims=1)
    Y_matrix = stack(solution.Y, dims=1)

    # Setup constraint check
    setupconstraint = (X_matrix .> MyAbsTol) .== Y_matrix
    if minimum(setupconstraint) == false
        rv -= 1
        push!(errmsg, "Invalid setup constraint!")
    end
    
    # Inventory computation
    inv = zeros(solution.data.P, solution.data.T)
    inv[:, 1] = solution.data.I0 + X_matrix[:, 1] - solution.data.D[:, 1]
    for i in 2:solution.data.T
        inv[:, i] = inv[:, i-1] + X_matrix[:, i] - solution.data.D[:, i]
    end
    
    inventoryconstraint = abs.(inv - I_matrix) .<= MyAbsTol
    if minimum(inventoryconstraint) == false
        rv -= 300
        deviation = maximum(abs.(inv - I_matrix))
        println(abs.(inv - I_matrix))
        push!(errmsg, "Invalid inventory constraint! Max deviation: $deviation")
    end

    # Cost validation
    cost_old = solution.Tcost
    solution.Tcost = update_Tcostv2(solution)
    if cost_old > solution.Tcost * (1 + MyAbsTol) || cost_old < solution.Tcost * (1 - MyAbsTol)
        push!(errmsg, "Incorrect total cost (diff: $(cost_old - solution.Tcost)) - corrected")
    end

    return rv, errmsg, setupconstraint, inventoryconstraint 
end

# check if job schedule is valid
function checkJobSchedule(solution::CLSPJSSolution)
    earliestStart = zeros(solution.data.No)
    latestFinish = zeros(solution.data.No)
    duration = zeros(solution.data.No)
    errmsg = 0

    for o in 1:solution.data.No
        i = solution.data.io[o]
        t = solution.data.lo[o]
        duration[o] = solution.X[i][t]*solution.data.po[o]+solution.Y[i][t]*solution.data.so[o]
        if abs(duration[o] - solution.schedule.duration[o])>MyAbsTol
            errmsg +=100
            println("Mismatch in duration of operation ", o)
        end
    end

    for o in solution.data.topological_order
        if solution.data.last_operation[o]
            t = solution.data.lo[o]
            m = solution.data.operation_to_machine[o]
            if t > 1
                earliestStart[o] = solution.data.cumulativeC[m][t-1]
            end
            if t == 1
                earliestStart[o] = 0.0
            end
        else
            earliestStart[o] = 0.0
        end
        if solution.data.has_pred[o]
            # predecessors finish
            for pred in solution.data.predecessors[o]
                if pred > 0
                    predFinish = earliestStart[pred] + duration[pred]
                    if predFinish > earliestStart[o]
                        earliestStart[o] = predFinish
                    end
                end
            end
        end
    end

    for o in solution.data.topological_order_reversed
        latestFinish[o] = solution.data.cumulativeC[solution.data.operation_to_machine[o]][solution.data.lo[o]]
        if solution.data.has_succ[o]
            # successors start
            for succ in solution.data.successors[o]
                if succ > 0
                    succStart = latestFinish[succ] - duration[succ]
                    if succStart < latestFinish[o]
                        latestFinish[o] = succStart
                    end
                end
            end
        end
    end

    for o in 1:solution.data.No
        if abs(earliestStart[o] - solution.schedule.earliestStart[o])>MyAbsTol
            errmsg +=20
            println("Mismatch in earliest start of operation ", o)
        end
        if abs(latestFinish[o] - solution.schedule.latestFinish[o])>MyAbsTol
            errmsg += 3
            println("Mismatch in latest finish of operation ", o)
        end
    end

    return errmsg
end

# reset demand list used in postponement routine
function reset_postponementDL!(alg::_2SCH)
    alg.postponementDL.amounts .= 0.0
    alg.postponementDL.itemIndex .= 0
    alg.postponementDL.periodIndex .= 0
    alg.nextPos = 1
end

# recalculate inventory cost
function update_Icost(solution::CLSPJSSolution)::Float64
    solution.Icost = 0
    for i in 1:solution.data.P
        sumInvP = 0
        for p in 1:solution.data.T
            sumInvP += solution.I[i][p]
        end
        solution.Icost += sumInvP * solution.data.h[i]
    end
    return solution.Icost
end

# recalculate setup cost
function update_Scost(solution::CLSPJSSolution)::Float64
    solution.Scost = 0
    for i in 1:solution.data.P
        sumSetupP = 0
        for p in 1:solution.data.T
            sumSetupP += solution.Y[i][p]
        end
        solution.Scost += sumSetupP * solution.data.s[i]
    end
    return solution.Scost
end

# recalculate production cost
function update_Pcost(solution::CLSPJSSolution)::Float64
    solution.Pcost = 0
    for i in 1:solution.data.P
        sumProdP = 0
        for p in 1:solution.data.T
            sumProdP += solution.X[i][p]
        end
        solution.Pcost += sumProdP * solution.data.pc[i]
    end
    return solution.Pcost
end

# recalculate tardiness penalty
function update_Ocost(solution::CLSPJSSolution)::Float64
    solution.Ocost = 0
    for operation in 1:solution.data.No
        tardiness = solution.schedule.earliestStart[operation] + solution.schedule.duration[operation] - solution.data.cumulativeC[solution.data.operation_to_machine[operation]][solution.data.lo[operation]]
        if tardiness > 0.0
            solution.Ocost += tardiness
        end
    end
    solution.Ocost *= 10000
    return solution.Ocost
end

# recalculate total cost
function update_Tcost(solution::CLSPJSSolution)::Float64
    solution.Tcost = update_Icost(solution) + update_Scost(solution) + update_Pcost(solution) + update_Ocost(solution)
    return solution.Tcost    
end 

# recalculate total cost in one loop (faster alternative for update_Tcost())
function update_Tcostv2(solution::CLSPJSSolution)::Float64
    solution.Icost = 0
    solution.Scost = 0
    solution.Pcost = 0
    solution.Ocost = 0
    @inbounds for i in 1:solution.data.P
        solution.Icost += @views sum(solution.I[i]) * solution.data.h[i]
        solution.Scost += @views sum(solution.Y[i]) * solution.data.s[i]
        solution.Pcost += @views sum(solution.X[i]) * solution.data.pc[i]
    end
    @inbounds for operation in 1:solution.data.No
        tardiness = solution.schedule.earliestStart[operation] + solution.schedule.duration[operation] - solution.data.cumulativeC[solution.data.operation_to_machine[operation]][solution.data.lo[operation]]
        if tardiness > 0.0
            solution.Ocost += tardiness
        end
    end
    solution.Ocost *= 10000
    solution.Tcost = solution.Icost + solution.Scost + solution.Pcost + solution.Ocost
    return solution.Tcost    
end

# try to add a lot (demand for period t) to a specific period (tt)
# returns: cost (incremental), maxamount (maximum amount that can be added)
function try_addLot2Period(solution::CLSPJSSolution, alg::_2SCH, i::Int64, t::Int64, tt::Int64, amount::Float64, dummyCost::Bool=false)::Tuple{Float64, Float64}
    activeShift = false # control if shift is enabled
    
    if amount < MyAbsTol
        maxamount = 0
        return 0, maxamount        
    end

    scost = 0   # setup cost (incremental)

    if solution.Y[i][tt] == false || alg.shiftProduction == true
        activeShift = true
    end

    # compute the amount that can be added
    maxamount = check_feasibility(solution, i, tt, amount)
    
    if !solution.Y[i][tt]      
        scost = solution.data.s[i]  # if there was no setup, there is additional setup cost
    end
    
    if maxamount <= MyAbsTol        # no production can be added to period tt
        maxamount = 0.0
        return -1, maxamount
    end
    
    #local variables to handle shifts
    cost_sav_inv = 0
    cost_sav_future = 0
    tt_future = tt + 1
    add_cap_inv = 0

    # ensure backup
    solution.schedule.duration_storage1 .= solution.schedule.duration
    solution.schedule.slack_storage1 .= solution.schedule.slack
    solution.schedule.earliestStart_storage1 .= solution.schedule.earliestStart
    solution.schedule.latestFinish_storage1 .= solution.schedule.latestFinish

    backupX = solution.X[i][tt] #backup the production amount
    backupY = solution.Y[i][tt]

    # simulation to update job schedule, must be reversed after shifts
    solution.X[i][tt] += maxamount
    solution.Y[i][tt] = true
    solution = change_duration!(solution, i, tt) # this solution has correct start and finish times, slack
    
    if activeShift
        # right shift
        if tt>1
            if (solution.I[i][tt-1] > MyAbsTol)
                maxfrominv = check_feasibility(solution, i, tt, solution.I[i][tt-1])
                add_cap_inv = maxfrominv   
                if maxfrominv > MyAbsTol
                    tt_back = tt-1
                    while maxfrominv > MyAbsTol
                        while solution.Y[i][tt_back] == false
                            cost_sav_inv += solution.data.h[i] * maxfrominv
                            tt_back -= 1
                        end
                        cost_sav_inv += solution.data.h[i] * maxfrominv
                        if maxfrominv - solution.X[i][tt_back] >= -MyAbsTol
                            maxfrominv -= solution.X[i][tt_back]
                            cost_sav_inv += solution.data.s[i]
                        else
                            break
                        end
                        tt_back -= 1
                    end
                end
            end
        end       
        solution.X[i][tt] += add_cap_inv
        solution = change_duration!(solution, i, tt)

        # left shift
        while tt_future <= solution.data.T
            if solution.Y[i][tt_future] == true
                amount_future = solution.X[i][tt_future]
                if dummyCost
                    amount_future += alg.cost_amounts[tt_future]
                end                
                avcap_future = check_feasibility(solution, i, tt, amount_future)
                if (avcap_future >= amount_future-MyAbsTol) && (amount_future*solution.data.h[i]*(tt_future-tt) < solution.data.s[i])
                    solution.X[i][tt] += amount_future
                    solution = change_duration!(solution, i, tt)
                    cost_sav_future += solution.data.s[i] - amount_future*solution.data.h[i]*(tt_future-tt)
                else
                    break
                end
            end
            tt_future += 1
        end
    end

    # reverse the changes
    solution.X[i][tt] = backupX
    solution.Y[i][tt] = backupY

    solution.schedule.duration, solution.schedule.duration_storage1 = solution.schedule.duration_storage1, solution.schedule.duration
    solution.schedule.slack, solution.schedule.slack_storage1 = solution.schedule.slack_storage1, solution.schedule.slack
    solution.schedule.earliestStart, solution.schedule.earliestStart_storage1 = solution.schedule.earliestStart_storage1, solution.schedule.earliestStart
    solution.schedule.latestFinish, solution.schedule.latestFinish_storage1 = solution.schedule.latestFinish_storage1, solution.schedule.latestFinish

    return scost + maxamount*solution.data.h[i]*(t-tt) - cost_sav_inv - cost_sav_future, maxamount
end

# add the amount (demand for period t) to production in period tt
# returns total cost of current plan
function addLot2Period(solution::CLSPJSSolution, alg::_2SCH, i::Int64, t::Int64, tt::Int64, amount::Float64)::CLSPJSSolution
    activeShift = false
        
    if amount <= MyAbsTol
        return solution
    end
        
    solution.X[i][tt] += amount
    solution.Pcost += amount * solution.data.pc[i]

    if solution.Y[i][tt] == false || alg.shiftProduction == true
        activeShift = true
    end
    if solution.Y[i][tt] == false
        solution.Y[i][tt] = true
        solution.Scost += solution.data.s[i]
    end

    solution = change_duration!(solution, i, tt)
        
    @inbounds for p in tt:t-1
        solution.I[i][p] += amount
    end
    solution.Icost += (t-tt) * solution.data.h[i] * amount

    solution.Ocost = 0.0
    @inbounds for operation in 1:length(solution.schedule.slack)
        tardiness = solution.schedule.earliestStart[operation] + solution.schedule.duration[operation] - solution.data.cumulativeC[solution.data.operation_to_machine[operation]][solution.data.lo[operation]]
        if tardiness > 0.0
            solution.Ocost += tardiness*10000
        end
    end

    solution.Tcost = solution.Icost + solution.Scost + solution.Pcost + solution.Ocost

    tt_future = tt + 1
    
    if activeShift
        # right shift
        if tt>1
            if (solution.I[i][tt-1] > MyAbsTol)
                maxfrominv = check_feasibility(solution, i, tt, solution.I[i][tt-1])
                if maxfrominv > MyAbsTol
                    solution.X[i][tt] += maxfrominv
                    solution = change_duration!(solution, i, tt)
                    tt_back = tt-1
                    while maxfrominv > MyAbsTol
                        while solution.Y[i][tt_back] == false
                            solution.I[i][tt_back] -= maxfrominv
                            solution.Icost -= solution.data.h[i] * maxfrominv
                            tt_back -= 1
                        end
                        solution.I[i][tt_back] -= maxfrominv
                        solution.Icost -= solution.data.h[i] * maxfrominv
                        if (maxfrominv - solution.X[i][tt_back])>= -MyAbsTol
                            maxfrominv -= solution.X[i][tt_back]
                            solution.Y[i][tt_back] = false
                            solution.Scost -= solution.data.s[i] 
                            solution.X[i][tt_back] = 0
                            solution = change_duration!(solution, i, tt_back)
                        else
                            solution.X[i][tt_back] -= maxfrominv
                            solution = change_duration!(solution, i, tt_back)
                            break
                        end
                        tt_back -= 1
                    end
                end
            end
        end
        
        # left shift
        while tt_future <= solution.data.T 
            if solution.Y[i][tt_future]==true
                avcap_future = check_feasibility(solution, i, tt, solution.X[i][tt_future])
                if (avcap_future >= solution.X[i][tt_future]-MyAbsTol && solution.X[i][tt_future]*solution.data.h[i]*(tt_future-tt) < solution.data.s[i])
                    solution.X[i][tt] += solution.X[i][tt_future]
                    solution = change_duration!(solution, i, tt)
                    solution.Y[i][tt_future] = false
                    solution.Scost -= solution.data.s[i]
                    @inbounds for ttt in tt:tt_future-1
                        solution.I[i][ttt] += solution.X[i][tt_future]
                    end
                    solution.Icost += solution.X[i][tt_future]*solution.data.h[i]*(tt_future-tt)   
                    solution.X[i][tt_future] = 0.0
                    solution = change_duration!(solution, i, tt_future)
                else
                    break
                end
            end
            tt_future += 1
        end
    end
    solution.Tcost = update_Tcostv2(solution)

    return solution
end

# case I of 2-SCH: use existing inventory to satisfy (partially) demand
# returns remaining amount to be scheduled
function _useInv(solution::CLSPJSSolution, dl::DemandList, alg::_2SCH, i::Int64, t::Int64, amount::Float64, pos::Int64)::Float64
    avCap = 0
    possible_amount = amount
    remaining_amount = 0

    if solution.I[i][t] > MyAbsTol     # if there is positive inventory
        @inbounds for p in t+1:solution.data.T
            av_amount = check_feasibility(solution, i, p, amount)
            if av_amount > MyAbsTol
                avCap += av_amount
            end
        end
        if avCap > MyAbsTol             # and if there is available capacity for production after period t
            if possible_amount > avCap
                possible_amount = avCap                               # compute what is the amount that can be satisfied from inventory: = d^new
            end
            if possible_amount > solution.I[i][t]
                possible_amount = solution.I[i][t]                    # compute what is the amount that can be satisfied from inventory: = d^new
            end

            remaining_amount = amount - possible_amount               # compute the remaining amount: =d^rem
            solution.I[i][t] -= possible_amount
            solution.Icost -= solution.data.h[i] * possible_amount
            solution.Tcost -= solution.data.h[i] * possible_amount
            alg.desactivatePostponement = true
            solution = addDemand(solution, dl, alg, i, t + 1, possible_amount, pos, true) # add the "new demand element" in the next period immediately
            if alg.stop && (!alg.replan)
                return 0.0
            else
                return remaining_amount
            end
        end
    end
    return amount
end

# options EXT and NEW: compute backwards cost of extending existing lots and/or creating new lots
# returns cost (incremental) of the cheapest option
function try_EXTandNEW(solution::CLSPJSSolution, alg::_2SCH, i::Int64, t::Int64, amount::Float64)::Float64
    NoSetupPeriodFound = false
    costEXT = 0
    costNEW = 1e+24
    thiscost = 1e+24
    tt = t-1
    ttt = t
    amount2 = 0
    possible_amount_rem = amount
    possible_amount_3 = 0
    alg.cost_amounts .= 0
    alg.postpone = 0

    # Backup before simulation
    solution.schedule.duration_storage3 .= solution.schedule.duration
    solution.schedule.slack_storage3 .= solution.schedule.slack
    solution.schedule.earliestStart_storage3 .= solution.schedule.earliestStart
    solution.schedule.latestFinish_storage3 .= solution.schedule.latestFinish     
    solution.X_storage3[i] .= solution.X[i]    #backup the production 


    while possible_amount_rem > MyAbsTol && tt >= 1           
        costEXT += possible_amount_rem * solution.data.h[i]

        amount2 = check_feasibility(solution, i, tt, possible_amount_rem)

        if solution.Y[i][tt] && amount2 > MyAbsTol
            possible_amount_rem -= amount2
            alg.cost_amounts[tt] = min(possible_amount_rem+amount2, amount2)
            # Simulate extension
            solution.X[i][tt] += amount2
            solution = change_duration!(solution, i, tt)
        end
        if ! solution.Y[i][tt]      
            # println("There is no setup in period ", tt) 
            if amount2 + MyAbsTol >= possible_amount_rem
                if ! NoSetupPeriodFound
                    thiscost, possible_amount_3 = try_addLot2Period(solution, alg, i, t, tt, possible_amount_rem, false)
                    if possible_amount_3 >= possible_amount_rem
                        if thiscost < costNEW
                            costNEW = thiscost
                            alg.cost_amounts[tt] = possible_amount_rem
                            ttt = tt
                            NoSetupPeriodFound = true
                        end
                    end
                end
            end
        end        
        tt -= 1
        if costEXT > costNEW
            break
        end
    end

    if possible_amount_rem > MyAbsTol
        costEXT = 1e+24
    end

    if alg.postponement && (costEXT != 1e+24)
        if (abs(costEXT-costNEW) <= alg.threshold * costEXT) && (alg.DL == 0) && ! alg.desactivatePostponement
            alg.postpone = 1
        end
    end

    # Undo simulation
    solution.X[i] .= solution.X_storage3[i]
    solution.schedule.duration, solution.schedule.duration_storage3 = solution.schedule.duration_storage3, solution.schedule.duration
    solution.schedule.slack, solution.schedule.slack_storage3 = solution.schedule.slack_storage3, solution.schedule.slack
    solution.schedule.earliestStart, solution.schedule.earliestStart_storage3 = solution.schedule.earliestStart_storage3, solution.schedule.earliestStart
    solution.schedule.latestFinish, solution.schedule.latestFinish_storage3 = solution.schedule.latestFinish_storage3, solution.schedule.latestFinish

    if costEXT <= costNEW
        alg.cost_amounts[ttt] = 0
        return costEXT
    else
        alg.cost_amounts[1:ttt-1] .= 0
        return costNEW
    end
end

# execute the cheapest option from EXT and NEW
function add_EXTorNEW(solution::CLSPJSSolution, alg::_2SCH, i::Int64, t::Int64)::CLSPJSSolution
    tt = t-1

    while tt >= 1
        if tt == 1
            solution = addLot2Period(solution, alg, i, t, tt, alg.cost_amounts[tt])
            return solution
        end
        solution = addLot2Period(solution, alg, i, t, tt, alg.cost_amounts[tt])
        tt -= 1
    end
end

# put demand element in the postponement demand list
function postponeDemand(solution::CLSPJSSolution, alg::_2SCH, i::Int64, t::Int64, amount::Float64, firstPostponement::Bool=true)::CLSPJSSolution
    if firstPostponement
        alg.postponementDL.itemIndex[alg.nextPos] = i
        alg.postponementDL.periodIndex[alg.nextPos] = t
        alg.postponementDL.amounts[alg.nextPos] = amount
        alg.nextPos += 1
    else
        for pos in 1:alg.nextPos-1
            if i == alg.postponementDL.itemIndex[pos] && t == alg.postponementDL.periodIndex[pos]
                alg.postponementDL.amounts[pos] = amount
            end
        end
        alg.keepInPostDL = true
    end
    alg.usePostDL = 0  
    return solution
end

# case IV of 2-SCH: use overtime (former), tardiness (here)
# returns total cost of current plan
function _useOvertime(solution::CLSPJSSolution, i::Int64, t::Int64, amount::Float64)::CLSPJSSolution
    tt = t-1
    t_inv = tt

    while tt >= 1 # fill available capacity with production in earlier periods as much as possible
        possible_amount = check_feasibility(solution, i, tt, amount)
        if possible_amount > MyAbsTol                                  
            amount2 = min(amount, possible_amount)
            if amount2 > MyAbsTol
                if solution.Y[i][tt] == false
                    solution.Y[i][tt] = true
                    solution.Scost += solution.data.s[i]
                end 
                solution.X[i][tt] += amount2
                solution = change_duration!(solution, i, tt)
                @inbounds for t_inv in tt:t-1
                    solution.I[i][t_inv] += amount2
                    solution.Icost += solution.data.h[i]*amount2*(t-t_inv)
                end
                amount -= amount2
            end
        end
        tt -= 1
    end
    if amount > MyAbsTol # put the rest in period t and accept tardiness with penalty
        if solution.Y[i][t] == false
            solution.Y[i][t] = true
            solution.Scost += solution.data.s[i]
        end   
        solution.X[i][t] += amount
        solution = change_duration!(solution, i, t)
        
        solution.Ocost = 0.0
        @inbounds for operation in 1:length(solution.schedule.slack)
            tardiness = solution.schedule.earliestStart[operation] + solution.schedule.duration[operation] - solution.data.cumulativeC[solution.data.operation_to_machine[operation]][solution.data.lo[operation]]
            if tardiness > 0.0
                solution.Ocost += tardiness*10000
            end
        end
    end

    solution.Tcost = update_Tcostv2(solution)
    return solution
end

# extend partial production plan: manages all cases of 2-SCH for a demand element
# returns total cost
function addDemand(solution::CLSPJSSolution, dl::DemandList, alg::_2SCH, i::Int64, t::Int64, amount::Float64, pos::Int64, activeUseInv::Bool=true, firstPostponement::Bool=true)::CLSPJSSolution    
    possible_amount = 0.0
    remaining_amount = 0.0
    cost = 0.0
    other_cost = 0.0

    alg.notAllNow = false

    if amount <= MyAbsTol   # check if amount is positive, otherwise no change of the production plan
        return solution
    end

    # case I: Use inventory
    if (activeUseInv && (t < solution.data.T))
        remaining_amount = _useInv(solution, dl, alg, i,t,amount, pos)
        if remaining_amount <= MyAbsTol
            return solution
        end
        alg.notAllNow = false
        amount = remaining_amount
    end

    # if we are in first period, we have to add the amount here
    if t == 1    
        solution = addLot2Period(solution, alg, i, t, t, amount)
        return solution
    end

    cost, possible_amount = try_addLot2Period(solution, alg, i,t,t,amount)    # compute option NOW
    remaining_amount = amount-possible_amount
    if possible_amount <= MyAbsTol
        alg.notAllNow = true
    end

    # case II: enough capacity, everything can be added to t
    if remaining_amount <= MyAbsTol  # if everything can be added but with additional cost (setup)
        other_cost = try_EXTandNEW(solution, alg, i, t, amount)
        if (cost <= other_cost) && ! alg.notAllNow  # option NOW - cheaper to add everything to t (either create a new lot or extend the existing lot)
            solution = addLot2Period(solution, alg, i, t, t, amount)
            return solution
        elseif (alg.DL == 0) && (alg.postpone == 1) # POSTPONEMENT routine - cheaper to create/extend previous lot/s but the difference between both is too small, so postpone 
            solution = postponeDemand(solution, alg, i, t, amount, firstPostponement)
            return solution
        else
            solution = add_EXTorNEW(solution, alg, i, t) # execute option EXT or NEW
            return solution
        end
    end

    # case III: not enough capacity, only part or nothing can be added to t
    # if there are no additional cost, but not the whole amount can be added
    if cost <= 0 && ! alg.notAllNow       
        solution = addLot2Period(solution,alg,i,t,t,possible_amount)    # add what is possible
        other_cost = try_EXTandNEW(solution,alg,i, t, remaining_amount) # try options EXT and NEW for the remaining amount (translates to options NOW+EXT and NOW+NEW options in case III)
        if other_cost == 1e+24                     # if options EXT and NEW are infeasible (=1e+24)
            solution = _useOvertime(solution, i, t, remaining_amount) # case IV: use overtime
            return solution
        end
        if (alg.DL == 0) && (alg.postpone == 1)   # POSTPONEMENT routine - options EXT and NEW are feasible but too small difference, so postpone
            solution = postponeDemand(solution, alg, i, t, remaining_amount, firstPostponement)
            return solution
        else  
            solution = add_EXTorNEW(solution, alg, i, t) # execute option EXT or NEW                                   
            return solution
        end

    # if there is additional cost and not everything can be added to t
    else
        # backup before simulation
        solution.schedule.duration_storage2 .= solution.schedule.duration
        solution.schedule.slack_storage2 .= solution.schedule.slack
        solution.schedule.earliestStart_storage2 .= solution.schedule.earliestStart
        solution.schedule.latestFinish_storage2 .= solution.schedule.latestFinish     
        solution.X_storage[i] .= solution.X[i]            #backup the production 
        solution.Y_storage[i] .= solution.Y[i]            #backup the setup 
        solution.I_storage[i] .= solution.I[i]            #backup the inventory

        solution = addLot2Period(solution, alg, i, t, t, possible_amount)  # simulation of option NOW to ensure correct job schedule                                          
        other_cost2 = try_EXTandNEW(solution, alg, i, t, remaining_amount) # try options EXT and NEW for the remaining amount (translates to options NOW+EXT and NOW+NEW options in case III)

        # undo simulation
        solution.X[i] .= solution.X_storage[i]
        solution.Y[i] .= solution.Y_storage[i]
        solution.I[i] .= solution.I_storage[i]
        solution.schedule.duration, solution.schedule.duration_storage2 = solution.schedule.duration_storage2, solution.schedule.duration
        solution.schedule.slack, solution.schedule.slack_storage2 = solution.schedule.slack_storage2, solution.schedule.slack
        solution.schedule.earliestStart, solution.schedule.earliestStart_storage2 = solution.schedule.earliestStart_storage2, solution.schedule.earliestStart
        solution.schedule.latestFinish, solution.schedule.latestFinish_storage2 = solution.schedule.latestFinish_storage2, solution.schedule.latestFinish
        solution.Tcost = update_Tcostv2(solution)

        other_cost = try_EXTandNEW(solution, alg, i, t, amount)  # try options EXT and NEW for the whole amount
        
        if (other_cost < other_cost2+cost) && (other_cost2 != 1e+24) # EXT or NEW is cheaper
            other_cost = try_EXTandNEW(solution, alg, i, t, amount)
        end
        if other_cost2 == 1e+24                              # if options EXT and NEW are infeasible (=1e+24)
            solution = _useOvertime(solution, i, t, amount)  # case IV: tardiness
            return solution
        end
        if (cost+other_cost2 <= other_cost) && ! alg.notAllNow # NOW+EXT or NOW+NEW is cheaper
            solution = addLot2Period(solution,alg,i,t,t,possible_amount)   # add what is possible
            other_cost2 = try_EXTandNEW(solution,alg,i, t, remaining_amount)
            if (alg.DL == 0) && (alg.postpone == 1)           # POSTPONEMENT routine - options EXT and NEW are feasible but too small difference, so postpone
                solution = postponeDemand(solution, alg, i, t, remaining_amount, firstPostponement)
                return solution
            else
                solution = add_EXTorNEW(solution, alg, i, t)      # execute option EXT or NEW
                return solution
            end
        elseif (alg.DL == 0) && (alg.postpone == 1)
            solution = postponeDemand(solution, alg, i, t, amount, firstPostponement)
            return solution
        else
            solution = add_EXTorNEW(solution, alg, i, t)          # execute option EXT or NEW
            return solution
        end
    end
end

# sort items according to some criteria
# returns a sorted list of item indices
function sortItems(solution::CLSPJSSolution, criterion::String)::Vector{Int64}
    @inbounds for i in 1:solution.data.P
        solution.itemlist[i] = i
    end
    if criterion == "SH"
        sort!(solution.itemlist, by= x -> (solution.data.s[x]/solution.data.h[x]),rev=true)
    elseif criterion == "HS"
        sort!(solution.itemlist, by = x -> (solution.data.h[x]/solution.data.s[x]), rev=true)
    elseif criterion == "TBO"
        sort!(solution.itemlist, by = x -> (solution.data.TBO[x]), rev=true)
    elseif criterion == "HSC"
        sort!(solution.itemlist, by = x -> (solution.data.h[x]/(solution.data.s[x]*solution.data.averageDemand[x])), rev=true)
    elseif criterion == "EC"
        sort!(solution.itemlist, by = x -> ((solution.data.s[x]/solution.data.TBO[x])+(solution.data.TBO[x]*solution.data.h[x]*solution.data.averageDemand[x]/2)), rev=true)
    elseif criterion == "ES"
        sort!(solution.itemlist, by = x -> ((solution.data.TBO[x]-1)*solution.data.s[x]-(solution.data.TBO[x]*(solution.data.TBO[x]-1)*solution.data.h[x]*solution.data.averageDemand[x]/2)), rev=true)
    elseif criterion == "ESC"
        sort!(solution.itemlist, by = x -> ((solution.data.TBO[x]-1)*solution.data.s[x]-(solution.data.TBO[x]*(solution.data.TBO[x]-1)*solution.data.h[x]*solution.data.averageDemand[x]/2))/solution.data.averageDemand[x], rev=true)
    elseif criterion == "CUI"
        sort!(solution.itemlist, by = x -> (sum(solution.data.D_Cap, dims=2)[x]), rev=true)
    elseif criterion == "S"
        sort!(solution.itemlist, by = x -> solution.data.s[x], rev=true)
    elseif criterion == "H"
        sort!(solution.itemlist, by = x -> solution.data.h[x], rev=true)
    elseif criterion == "INCI"
        return solution.itemlist
    else
        print("Invalid Criterion!")
    end
    return solution.itemlist
end

# sort items according to some criteria
# returns a sorted list of period indices
function sortPeriods(solution::CLSPJSSolution, criterion::String)::Vector{Int64}
    @inbounds for p in 2:solution.data.T
        solution.periodlist[p-1] = p
    end
    if criterion == "INCP"
        return solution.periodlist
    elseif criterion == "DECP"
        sort!(solution.periodlist, by = x -> x, rev=true)
    elseif criterion == "CUP"
        sort!(solution.periodlist, by = x -> (sum(solution.data.D_Cap, dims=1)[x]), rev=true)
    else
        print("Invalid Criterion!")
    end
    return solution.periodlist
end

# generate demand list with specific sorting critieria (only for cost-based criteria)
# returns demand list
function generateDL(dl::DemandList, solution::CLSPJSSolution, criterionPeriod::String, criterionItem::String, choice::String)::DemandList
    solution.itemlist .= sortItems(solution, criterionItem)                  
    solution.periodlist .= sortPeriods(solution, criterionPeriod)           
    
    # reset demand list
    dl.itemIndex .= 0
    dl.periodIndex .= 0
    dl.amounts .= 0

    pos = 1
    if choice == "PI"
        @inbounds for t in solution.periodlist, i in solution.itemlist
            D_it = solution.data.D[i, t]
            if D_it > MyAbsTol
                dl.itemIndex[pos] = i
                dl.periodIndex[pos] = t
                dl.amounts[pos] = D_it
                pos += 1
            end
        end
    elseif choice == "IP"
        @inbounds for i in solution.itemlist, t in solution.periodlist
            D_it = solution.data.D[i, t]
            if D_it > MyAbsTol
                dl.itemIndex[pos] = i
                dl.periodIndex[pos] = t
                dl.amounts[pos] = D_it
                pos += 1
            end
        end
    else
        print("Invalid choice!")
    end
    return dl
end

# generate demand list with specific sorting critieria (only for capacity-based criteria)
# returns demand list
function generateDLdemand(dl::DemandList, solution::CLSPJSSolution, criterionPeriod::String, criterionDemand::String)::DemandList
    # reset demand list
    dl.itemIndex .= 0
    dl.periodIndex .= 0
    dl.amounts .= 0

    solution.periodlist = sortPeriods(solution, criterionPeriod)
    pos = 1
    @inbounds for p in solution.periodlist
        @inbounds for i in 1:solution.data.P
            solution.itemlist[i] = i
        end 
        if criterionDemand == "D"
            sort!(solution.itemlist, by = x -> solution.data.D_Cap[x, p], rev=true)
        elseif criterionDemand == "DH"
            sort!(solution.itemlist, by = x -> solution.data.D_Cap[x, p] / solution.data.h[x])
        elseif criterionDemand == "DS"
            sort!(solution.itemlist, by = x -> solution.data.D_Cap[x, p] / solution.data.s[x])
        elseif criterionDemand == "DHS"
            sort!(solution.itemlist, by = x -> solution.data.D_Cap[x, p] * solution.data.h[x] / solution.data.s[x])
        elseif criterionDemand == "DSH"
            sort!(solution.itemlist, by = x -> solution.data.D_Cap[x, p] / (solution.data.h[x] * solution.data.s[x]))
        else
            print("Invalid Criterion!")
        end
        @inbounds for i in solution.itemlist
            D_it = solution.data.D[i, p]
            if D_it > MyAbsTol
                dl.itemIndex[pos] = i
                dl.periodIndex[pos] = p
                dl.amounts[pos] = D_it
                pos += 1
            end
        end
    end
    return dl
end

# plan a demand list from pos1 till pos2
# returns total cost of production plan
function planDemandList(solution::CLSPJSSolution, dl::DemandList, alg::_2SCH, resetSol::Bool = true, pos1::Int64 = 1, pos2::Int64 = 1000000000, resetPostDL::Bool = true)::Tuple{Float64, CLSPJSSolution}
    pos2 = min(pos2, solution.data.nonzerodemand)
    if resetSol                # reset production plan if true, otherwise continue with plan
        solution = resetSolution!(solution)
    end
    if resetPostDL
        reset_postponementDL!(alg)
    end

    @inbounds for pos in pos1:pos2 # go through main demand list
        if alg.stop && (!alg.replan)
            reset_postponementDL!(alg)
            break
        end
        alg.usePostDL = 1
        alg.desactivatePostponement = false
        solution = addDemand(solution, dl, alg, dl.itemIndex[pos], dl.periodIndex[pos], dl.amounts[pos], pos) # add demand element to the partial production plan
        if solution.Ocost > MyAbsTol
            alg.stop = false # ->>>>>> if false, do not replan with different holding cost, i.e., no cost adaptation
        end

        if alg.stop && (!alg.replan)
            reset_postponementDL!(alg)
            break
        end

        @inbounds for pos_in_DL2 in 1:alg.nextPos-1 # go through postponement demand list
            alg.keepInPostDL = false
            if alg.usePostDL && (dl.itemIndex[pos] == alg.postponementDL.itemIndex[pos_in_DL2])
    
                solution = addDemand(solution, dl, alg, alg.postponementDL.itemIndex[pos_in_DL2], alg.postponementDL.periodIndex[pos_in_DL2], alg.postponementDL.amounts[pos_in_DL2], pos, true, false)
                if solution.Ocost > MyAbsTol
                    alg.stop = false # ->>>>>> if false, do not replan with different holding cost, i.e., no cost adaptation
                end
                if ! alg.keepInPostDL
                    alg.postponementDL.itemIndex[pos_in_DL2] = -1 # "delete" this element
                end

                if alg.stop && (!alg.replan)
                    reset_postponementDL!(alg)
                    break
                end
            end
        end
    end

    # after we went through the whole main demand list, ensure that nothing remained unscheduled in the postponement demand list
    alg.DL = 1 
    @inbounds for posDL2 in 1:alg.nextPos-1
        if alg.postponementDL.itemIndex[posDL2] != -1
            solution = addDemand(solution, dl, alg, alg.postponementDL.itemIndex[posDL2], alg.postponementDL.periodIndex[posDL2], alg.postponementDL.amounts[posDL2], solution.data.nonzerodemand, true, false)
            if solution.Ocost > MyAbsTol
                alg.stop = false # ->>>>>> if false, do not replan with different holding cost, i.e., no cost adaptation
            end
        end
    end
    alg.DL = 0 # reset indicator which DL is used
    return solution.Tcost, solution
end

# write solution to a text file
function CLSP_write_complete_solution(sol::CLSPJSSolution, path::String, rule::String, name::String="")
    if name == ""
        filename = joinpath(path, string(sol.data.name,"_",rule,"_fullsol.txt"))
    else
        filename = joinpath(path, string(name,"_",rule,"_fullsol.txt"))
    end

    open(filename,"w") do file
        println(file,"Name:\t",sol.data.name)
        println(file,"Items:\t",sol.data.P)
        println(file,"Periods:\t",sol.data.T)
        print(file, "HoldingCost:\t")
        writedlm(file,sol.data.h')
        print(file, "SetupCost:\t")
        writedlm(file,sol.data.s')
        println(file,"Demand:")
        writedlm(file,sol.data.D)
        print(file, "InitialInventory:\t")
        writedlm(file,sol.data.I0')
        print(file, "Capacity:\t")
        writedlm(file,sol.data.C')
        print(file, "ProductionCoefficient:\t")
        writedlm(file,sol.data.a')
        print(file, "SetupCoefficient:\t")
        writedlm(file,sol.data.b')
        println(file,"TotalCost:\t",sol.Tcost)
        println(file,"ProductionCost:\t",sol.Pcost)
        println(file,"InventoryCost:\t",sol.Icost)
        println(file,"SetupCost:\t",sol.Scost)
        println(file,"TardinessPenalty:\t",sol.Ocost)
        println(file,"Setup:")
        writedlm(file,sol.Y)
        println(file,"Lotsizes:")
        writedlm(file,sol.X)
        println(file,"Inventory:")
        writedlm(file,sol.I)
        println(file,"Duration:\t")
        writedlm(file,sol.schedule.duration)
        println(file,"Earliest Start Time:\t")
        writedlm(file,sol.schedule.earliestStart)
        println(file,"Latest Finish Time:\t")
        writedlm(file,sol.schedule.latestFinish)
        println(file,"Slack:\t")
        writedlm(file,sol.schedule.slack)
        println(file,"Predecessors:\t")
        writedlm(file,sol.data.predecessors)
        println(file,"Successors:\t")
        writedlm(file,sol.data.successors)
    end
end                 

# run 2SCH for all Wolosewich instances and save output
function _run2SCH_WolosewichInstances()
    path = "C:/Users/ddziuba/Documents/CLSP with JS Scheduling/Instances Wolosewich/"
    job_shop_size = ["6x6", "10x10", "20x5"]
    T_values = [2, 3, 5, 10, 20, 30, 50, 100]

    variations = ["Defaut", "Variation_capacite", "Variation_cout_setup", "Variation_demande", "Variation_lead_time"]
    # "Defaut" -> 1, "Variation_capacite" -> 4, "Variation_cout_setup" -> 4, "Variation_demande" -> 5, "Variation_lead_time" -> 4

    # Output Dataframe
    output = DataFrame(Instance = String[], Rule = String[], CPUTime = Float64[], Tcost = Float64[], Pcost = Float64[],
                        Icost = Float64[], Scost = Float64[], Ocost = Float64[],
                        errmsg = Int64[], errmsgText = Vector[], schedule_errmsg= Int64[])

    # Criteria for sorting rules
    cost_crit = ["CUI", "EC","ES","ESC", "H", "HS", "HSC", "S","SH","TBO"]
    period_crit = ["CUP", "DECP","INCP"]
    capacity_crit = ["D","DH", "DHS", "DS", "DSH"]
    key_crit = ["IP", "PI"]

    # All sorting rules (75)
    sorting_rules = []
    for pc in period_crit, ic in cost_crit, k in key_crit
        rule = pc*"-"*ic*"-"*k
        push!(sorting_rules, rule)
    end
    for pc in period_crit, dc in capacity_crit
        rule = pc*"-"*dc
        push!(sorting_rules, rule)
    end

    # Run 2SCH for all instances
    for size in job_shop_size[1:end], T in T_values[1:end], var in variations[1:end], i in 1:5
        instance_path = path*size*"/"*"T="*string(T)*"/"*var*"/"*string(i)*"/"
        if isdir(instance_path)
            filename = size*"-"*"T="*string(T)*"-"*var*"-"*string(i)
            data = W_read_CLSPJSData(instance_path)                       # read instance
            for pc in period_crit, ic in cost_crit, k in key_crit
                alg = _2SCH(data)                                         # set up 2-SCH version
                solution =CLSPJSSolution(data)                            # set up solution container
                dl = DemandList(data.nonzerodemand)                       # set up demand list
                rule = pc*"-"*ic*"-"*k
                time = @elapsed begin                                      # time the following part
                    dl = generateDL(dl, solution, pc, ic, k)               # Step 1 - generate demand list
                    obj, solution = planDemandList(solution, dl, alg)      # Step 2 - plan demand list           
                end
                push!(output, [instance_path, rule, time, obj, solution.Pcost, solution.Icost, solution.Scost, solution.Ocost,
                                        checkSolution(solution)[1], checkSolution(solution)[2], checkJobSchedule(solution)])

                alg = nothing
                solution = nothing
                dl = nothing
            end

            for pc in period_crit, dc in capacity_crit
                alg = _2SCH(data)                                         # set up 2-SCH version
                solution =CLSPJSSolution(data)                            # set up solution container
                dl = DemandList(data.nonzerodemand)                       # set up demand list
                rule = pc*"-"*dc
                time = @elapsed begin                                   # time the following part
                    dl = generateDLdemand(dl, solution, pc, dc)         # Step 1 - generate demand list
                    obj, solution = planDemandList(solution, dl, alg)   # Step 2 - plan demand list               
                end
                push!(output, [instance_path, rule, time, obj, solution.Pcost, solution.Icost, solution.Scost, solution.Ocost, 
                        checkSolution(solution)[1], checkSolution(solution)[2], checkJobSchedule(solution)])
                
                alg = nothing
                solution = nothing
                dl = nothing
            end
            data = nothing
        end
    end
    CSV.write("output.csv", output)
end

# run 2SCH for one instance from Wolosewich et al. and save output
function run2SCH_oneInstance(instance_path::String)
    # Output Dataframe
    output = DataFrame(Instance = String[], Rule = String[], CPUTime = Float64[], Tcost = Float64[], Pcost = Float64[],
                        Icost = Float64[], Scost = Float64[], Ocost = Float64[],
                        errmsg = Int64[], errmsgText = Vector[], schedule_errmsg= Int64[])

    # Criteria for sorting rules
    cost_crit = ["CUI", "EC","ES","ESC", "H", "HS", "HSC", "S","SH","TBO"]
    period_crit = ["CUP", "DECP","INCP"]
    capacity_crit = ["D","DH", "DHS", "DS", "DSH"]
    key_crit = ["IP", "PI"]

    # All sorting rules (75)
    sorting_rules = []
    for pc in period_crit, ic in cost_crit, k in key_crit
        rule = pc*"-"*ic*"-"*k
        push!(sorting_rules, rule)
    end
    for pc in period_crit, dc in capacity_crit
        rule = pc*"-"*dc
        push!(sorting_rules, rule)
    end

    # Run 2SCH for instance at path
    if isdir(instance_path)
        data = W_read_CLSPJSData(instance_path)                       # read instance
        for pc in period_crit, ic in cost_crit, k in key_crit
            alg = _2SCH(data)                                         # set up 2-SCH version
            solution =CLSPJSSolution(data)                            # set up solution container
            dl = DemandList(data.nonzerodemand)                       # set up demand list
            rule = pc*"-"*ic*"-"*k
            time = @elapsed begin                                      # time the following part
                dl = generateDL(dl, solution, pc, ic, k)               # Step 1 - generate demand list
                obj, solution = planDemandList(solution, dl, alg)      # Step 2 - plan demand list           
            end
            push!(output, [instance_path, rule, time, obj, solution.Pcost, solution.Icost, solution.Scost, solution.Ocost,
                                    checkSolution(solution)[1], checkSolution(solution)[2], checkJobSchedule(solution)])

            alg = nothing
            solution = nothing
            dl = nothing
        end

        for pc in period_crit, dc in capacity_crit
            alg = _2SCH(data)                                         # set up 2-SCH version
            solution =CLSPJSSolution(data)                            # set up solution container
            dl = DemandList(data.nonzerodemand)                       # set up demand list
            rule = pc*"-"*dc
            time = @elapsed begin                                   # time the following part
                dl = generateDLdemand(dl, solution, pc, dc)         # Step 1 - generate demand list
                obj, solution = planDemandList(solution, dl, alg)   # Step 2 - plan demand list               
            end
            push!(output, [instance_path, rule, time, obj, solution.Pcost, solution.Icost, solution.Scost, solution.Ocost, 
                    checkSolution(solution)[1], checkSolution(solution)[2], checkJobSchedule(solution)])
            
            alg = nothing
            solution = nothing
            dl = nothing
        end
        data = nothing
    end
    CSV.write("output.csv", output)
end

# _run2SCH_WolosewichInstances()

# run2SCH_oneInstance("C:/Users/ddziuba/Documents/CLSP with JS Scheduling/Instances Wolosewich/10x10/T=10/Defaut/1/")

end