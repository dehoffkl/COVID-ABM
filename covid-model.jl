using Agents
using BenchmarkTools
using Distributions
using Dates
using Random

"""
    CovidAgent

Definition of agent structure - single person in model

Inherits from Agents.AbstractAgent

Contains the following elements:
    - id::Int: Identifier of agent
    - pos::NTuple{2, Float64}: Current position within 2D space
    - age::Int: Age of agent
    - status::Symbol: Either :S, :I, or :R
    - infect_t::Float64: Time since initial infection
    - severity::Float64: Current infection severity (out of 100)
    - crit_t::Float64: Amount of time with severity over model defined Scrit
    - recover_t::Float64: Time since agent is considered recovered
    - mask::Bool: Whether agent wears a mask
    - β::Float64: Current transmission probability
    - βparams::Dict: Parameters for updating β during infection
    - Sparams::Dict: Parameters for updating severity during infection
    - r::Float: reinfection probability
    - quarantined::Bool: Whether infection has been discovered and the agent quarantined
    - hospitalized::Bool: Whether agent has been put into hospital due to infection
    - ignore_quarantine::Bool: Whether agent is ignoring stay-at-home orders
    - vel::NTuple{2, Float64}: Agent's current velocity as a vector
    - mass::Float64: Agent's current mass (prevents movement during collisions)
"""
mutable struct CovidAgent <: AbstractAgent
    id::Int
    pos::NTuple{2,Float64}
    age::Int
    status::Symbol
    infect_t::Float64
    severity::Float64
    crit_t::Float64
    recover_t::Float64
    mask::Bool
    β::Float64
    βparams::Dict
    Sparams::Dict
    r::Float64
    quarantined::Bool
    hospitalized::Bool
    ignore_quarantine::Bool
    vel::NTuple{2,Float64}
    mass::Float64
end

"""
   gen_β(min_Kβ, max_kβ, βmax)

    Generates values for βparams for an agent.  Many of these
    values are based on observed data from the pandemic.  The
    returned dictionary contains the following parameters

     - tβmax::Float64: Time at max contagion, 2-6 days by observed cycle
     - tβmed::Float64: Time at median contagion, halfway between 2 days and tβmax
     - Kβ::Float64: How quickly the transmission probability builds (arbitrary between 5 and 10)
     - βmax::Float64: Maximum transmission probability (Normal distribution around βmax with 0.05*βmax std dev)
     - η::Float64: How quickly transmission probability decreases (not contagious by 11 days)

    # Arguments
     - min_Kβ:Float64: Minimum value for Kβ
     - max_Kβ:Float64: Maximum value for Kβ
     - βmax:Float64: Average transmission probability (in percentage)
"""
function gen_β(Sparams, min_Kβ=5, max_kβ=10, βmax=75)
    lb = min_Kβ < max_kβ ? min_Kβ : max_kβ;
    ub = min_Kβ > max_kβ ? min_Kβ : max_kβ;
    ub = min_Kβ == max_kβ ? max_kβ : max_kβ + 1;
    βparams = Dict() 
    #Looking for a tight distribution around 5% (± 2%)
    test = Sparams["tS0"] - rand(Uniform(1., 3.))
    βparams["tβmax"] = test > 0.0 ? test : 0.0
    βparams["tβmed"] = 2. + (βparams["tβmax"] - 2.)/2.
    βparams["Kβ"] = rand(Uniform(lb, ub))
    βparams["βmax"] = rand(Normal(βmax, (βmax*0.05)))/100
    βparams["η"] = βparams["βmax"] / (11-βparams["tβmax"])
    return βparams
end

"""
   gen_S(age)

    Generates values for Sparams for an agent.  Most of these values are arbitrary

    - tSmax::Float64: Time at max severity, 1-14 days after onset of symptoms
    - tS0::Float64: Time at onset of symptoms, 2-14 days after infection
    - KS::Float64: How quickly the severity builds (arbitrary between 0.1 and 2)
    - Smax::Float64: Maximum severity (limited normal around logistic function)
    - γ::Float64: How quickly severity decreases (linear recovery)
"""
function gen_S(age)
    Sparams = Dict()
    #looking for a loose distribution around logistic function
    #assume a normal distribution with variance 20.  Values
    #outside [0, 100] will be truncated to the range limit
    s_max_mean = 100/(1+exp(-0.2*(age-70)))
    tmp = rand(Normal(s_max_mean, 20))
    tmp > 100 ? Sparams["Smax"] = 100 : 
        tmp < 0 ? Sparams["Smax"] = 0 :
            Sparams["Smax"] = tmp
    Sparams["tS0"] = rand(Uniform(2, 14))
    Sparams["tSmax"] = Sparams["tS0"] + rand(Uniform(1, 14))
    Sparams["Ks"] = rand(Uniform(0.1, 2))
    Sparams["γ"] = rand(Normal(10, 2))
    return Sparams
end

function dayFrac(p)
    return convert(Float64, convert(Dates.Millisecond, p) / convert(Dates.Millisecond, Day(1)))
end

"""
    initialize([extents, speed, mask_perc, ageRange, percInfected,
         percIgnore, dt, initial_qlevel, interaction_radius,
         detection_threshold, critical_threshold)

    Initialize the ABM.  Creates the model space and the agents

    # Arguments
    - extents:NTuple{2, Float64}: Dimensions of the model space
    - N:Int64: Number of agent in the model
    - speed:Float64: Absolute speed of agents
    - mask_perc:Float64: Percentage of agents that use mask
    - ageRange:NTuple{2, Int64}: Range of ages for agents
    - percInfected:Float64: Initial percentage of infected agents
    - percIgnore:Float64: Percentage of agents ignoring stay-at-home orders
    - dt:Float64: Model time step (in days)
    - initial_qlevel:Int64: Initial stay-at-home order level (0 = maximum, 4 = fully open)
    - interaction_radius:Float64: How close two agents need to be in order to potentially transmit
    - detection_threshold:Float64: Threshold for detecting infection
    - critical_threshold:Float64: Threshold for being considered seriously ill
    - reinfect_prob:Float64: Maximum probability for reinfection
    - init_date:Date: Initial date for model
"""
function initialize(; extents=(10,10),
    N = 500,
    speed=0.002, 
    mask_perc=0.9,
    ageRange = (10, 90),
    percInfected = 0.1,
    percIgnore = 0.2,
    dt = Dates.Day(1),
    initial_qlevel = 0,
    interaction_radius = 2,
    detection_threshold = 10,
    critical_threshold = 75,
    reinfect_prob = 0.05,
    init_date = Date(2020, 03, 01),
    Kβ_min = 5,
    Kβ_max = 10,
    βmax = 75,
    mask_effect = 0.67,
    rseed=-1
    )

    if rseed > 0
        Random.seed!(rseed)
    end
    properties = Dict(
        :dt => dayFrac(dt),
        :qlevel => initial_qlevel,
        :radius => interaction_radius,
        :Sdetect => detection_threshold,
        :Scrit => critical_threshold,
        :speed => speed,
        :reinf => reinfect_prob,
        :tick => convert(DateTime, init_date),
        :Kβ_min => Kβ_min,
        :Kβ_max => Kβ_max,
        :βmax => βmax,
        :mask_effect => mask_effect,
        :tick_delt => dt,
    )

    space2d = ContinuousSpace(2; periodic=true, extend=extents)
    model = ABM(CovidAgent, space2d, 
                properties = properties)
    ageDist = DiscreteUniform(ageRange...)


    #add agents to the model
    for ind in 1:N
        age = rand(ageDist)
        mask = rand() > mask_perc
        sym = rand() < percInfected ? :I : :S
        ignore = rand() < percIgnore
        mass = ignore ? 1.0 : Inf
        vel = ignore ? calc_speed(speed) : (0.0, 0.0)
        Sparams = gen_S(age)
        βparams = gen_β(Sparams, Kβ_min, Kβ_max, βmax)
        add_agent!(model, 
            age, #agent age
            sym, #current infection status
            0, #infection elapsed time
            0, #infection severity
            0, #time over critical threshold
            0, #time since entered recovery status
            mask, #whether mask is used
            0, #current transmission probability
            βparams, #severity parameters
            Sparams, #transmission probability params
            0.01, #reinfection probability
            false, #quarantined
            false, #hospitalized
            ignore, #ignoring quarantine
            vel, #current velocities
            mass) #current mass
    end
    index!(model)
    return model
end

"""
    agent_step!(agent, model)

Update a single agent for a single time step.

This function will first move the agent in space, then update
the transmission probability β and severity.  It will then
check to see if the agent is infected and quarantined, hospitalized, or passed.
"""
function agent_step!(agent, model)
    move_agent!(agent, model, model.dt)
    update_agent!(agent, model)
    check_status!(agent, model)
end

"""
    update_agent!(a, m)

This function will update an agent.  If the agent is infected, it will
calculate the transmission probability and severity at the new time step.
If the severity drops below 0, then the agent is considered recovered (:R)
immediately
"""
function update_agent!(a, m)
    if a.status == :I
        a.infect_t += m.dt

        if a.severity > m.Scrit
            a.crit_t += m.dt
        end

        t = a.infect_t
        #update β
        βp = a.βparams
        Δβ = 0
        if t < βp["tβmax"]
            Δβ = βp["βmax"]*βp["Kβ"]*(exp(-βp["Kβ"]*(t-βp["tβmed"]))) / 
                (1+exp(-βp["Kβ"]*(t-βp["tβmed"])))^2
        else
            Δβ = -βp["η"]
        end
        a.β += Δβ*m.dt
        if a.β < 0 
            a.β = 0
        end

        #update severity
        Sp = a.Sparams
        ΔS = 0
        if t < Sp["tSmax"]
            ΔS = Sp["Smax"]*Sp["Ks"]*(exp(-Sp["Ks"]*(t-Sp["tS0"]))) / 
                (1+exp(-Sp["Ks"]*(t-Sp["tS0"])))^2
            if a.hospitalized 
                ΔS = 0.5*ΔS
            end
            #print(ΔS)
        else
            ΔS = -Sp["γ"]
        end
        a.severity += ΔS*m.dt
        if a.severity <= 0
            a.status = :R
            a.recover_t = 0
            reset_agent!(a,m)
        end
    elseif a.status == :R
        a.recover_t += m.dt
        if a.r < m.reinf
            a.r += m.dt * rand(Normal(m.reinf, m.reinf/6.))
        end
    end
end

"""
    check_status(a, m)

This function will check the status of an infected agent.  If
the infection severity is above the model detection threshold, 
or a random number  is less than the model severity, then the 
agent is quarantined.  We then check to see if the agent should be
hospitalized.  If the severity is greater than the critical threshold
and a random number is less than 1/4 the current severity, then the
patient is sent to the hospital.  If the patient has had a severity greater
than the critial threshold for 5 days, then the patient is considered 
to have died.  If the time since initial infection is greater than 30
days, then the agent is considered to be recovered.
"""
function check_status!(a, m)
    if a.status == :I
        #check for detection and quarantine
        if a.severity > m.Sdetect || rand() < a.severity
            a.quarantined == true
            a.hospitalized = false
            a.mass = Inf
            a.vel = (0.0,0.0)
        end

        #check for hospitalization
        if a.severity > m.Scrit && rand() < a.severity / 4.
            a.quarantined = true
            a.hospitalized = true
            a.mass = Inf
            a.vel = (0.0,0.0)
        end

        #check for death
        if a.crit_t > 3
            kill_agent!(a, m)
        end

        #check for recovery
        if a.infect_t > 30
            a.status = :R
            a.recover_t = 0
            a.r = 0
            reset_agent!(a, m)
        end
    end
end

"""
    reset_agent!(a, m)

This function resets the agent after recovery with a new set of infection
parameters.  It releases the agent from quarantine and from the hospital.
"""
function reset_agent!(a, m)
    a.β = 0
    a.severity = 0
    a.infect_t = 0
    a.crit_t = 0
    a.quarantined = false
    a.hospitalized = false
    a.mass = a.ignore_quarantine ? 1.0 : Inf
    a.ignore_quarantine ? a.vel = calc_speed(m.speed) : a.vel = (0.0, 0.0)
    a.Sparams = gen_S(a.age)
    a.βparams = gen_β(a.Sparams, m.Kβ_min, m.Kβ_max, m.βmax)
end

"""
    transmit!(a1, a2, m)

This function checks for transmission between two agents.  It first checks that
only one of a1, a2 are infected.  Then it adjusts the β (probability of transmission)
based on mask usage by both parties.  If the healthy party has an :S status, then
it uses the transmission probability β from the infected individual.  If the
healty party has the :R status, then it uses the reinfection probability r from
the healthy individual 
"""
function transmit!(a1, a2, m)
    # for transmission, only 1 can have the disease (otherwise nothing happens)
    count(a.status == :I for a in (a1, a2)) ≠ 1 && return
    infected, healthy = a1.status == :I ? (a1, a2) : (a2, a1)

    actual_β = infected.β
    if healthy.status == :R
        actual_β = healthy.r
    end

    healthy.mask ? actual_β *= m.mask_effect : actual_β = actual_β
    infected.mask ? actual_β *= m.mask_effect : actual_β = actual_β

    rand() > actual_β && return

    healthy.infect_t = 0
    healthy.status = :I
end

"""
    calc_speed(max_speed, multipler)

Calculates the velocity of an agent based on a maximum speed defined in 
the model and a multiplier defined by the current quarantine level
"""
function calc_speed(max_speed, multiplier = 1.0)
    return sincos(2π * rand()) .* (multiplier * max_speed)
end

"""
    change_quarantine!(model)

This function changes the population quarantine level.  At level 0,
those who follow the quarantine do not move (isolated).  At level 1,
we increase the speed of those following quarantine with age < 70 
and are not specifically quarantined or hospitalized to 33% of model speed.  
At level 2, we repeat the speed increase, and set the appropriate agents
to 66% of model speed.  At level 3, we repeat again, and set the appropriate
agents to 100% of model speed.  At level 4, we release all agents that are
not specifically quarantined or hospitalized to full speed
"""
function change_quarantine!(model)
    if model.qlevel == 0
        for a in allagents(model)
            if !a.ignore_quarantine
                a.vel = (0,0)
                a.mass = Inf
            end
        end

    elseif model.qlevel == 1
        for a in allagents(model)
            if a.age < 70
                if !(a.quarantined || a.hospitalized) && !a.ignore_quarantine
                    a.vel = calc_speed(model.speed, 0.33)
                    a.mass = 1.0 
                end
            end
        end
    elseif model.qlevel == 2
        for a in allagents(model)
            if a.age < 70
                if !(a.quarantined || a.hospitalized) && !a.ignore_quarantine
                    a.vel = calc_speed(model.speed, 0.67)
                    a.mass = 1.0 
                end
            end
        end
    elseif model.qlevel == 3
        for a in allagents(model)
            if a.age < 70
                if !(a.quarantined || a.hospitalized) && !a.ignore_quarantine
                    a.vel = calc_speed(model.speed)
                    a.mass = 1.0 
                end
            end
        end
    elseif model.qlevel == 4
        for a in allagents(model)
            if !(a.quarantined || a.hospitalized) && !a.ignore_quarantine
                a.vel = calc_speed(model.speed)
                a.mass = 1.0 
            end
        end
    end
end


"""
    model_step!(model)

This runs a single step of the model.  It tests for interacting agents (within
interaction_radius units), and attempts to transmit the virus.  It then adjusts
speed and direction of interacting agents through elastic collisions.
"""
function model_step!(model)
    r = model.radius
    for (a1, a2) in interacting_pairs(model, r, :nearest)
        transmit!(a1, a2, model)
        elastic_collision!(a1, a2, :mass)
    end

    if model.tick == Dates.DateTime(2020, 3, 26)
        model.qlevel = 0
        change_quarantine!(model)
    elseif model.tick == Dates.DateTime(2020, 5, 4)
        model.qlevel = 1
        change_quarantine!(model)
    elseif model.tick == Dates.DateTime(2020, 6, 5)
        model.qlevel = 2
        change_quarantine!(model)
    elseif model.tick == Dates.DateTime(2020, 6, 19)
        model.qlevel = 3
        change_quarantine!(model)
    elseif model.tick == Dates.DateTime(2020, 7, 5)
        model.qlevel = 4
        change_quarantine!(model)
    end

    model.tick += model.tick_delt
end

