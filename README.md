# COVID-ABM

This is an agent-based model meant to simulate the spread of COVID-19 through a population.

## Major Model Parameters
|Data|Type|Description|Expected Values|
|----|----|----|----|
|N|Int|Number of agents|10K-100K|
|<img src="https://render.githubusercontent.com/render/math?math=S_{detect}">|Float|Detection severity threshold|5-20|
|<img src="https://render.githubusercontent.com/render/math?math=S_{crit}">|Float|Criticial severity threshold|60-80|
|(X,Y)|Tuple(Float, 2)|Dimensions of model area|(1-1000, 1-1000)|


## Major Agent Parameters
I will be creating a population where each member will have the following information:

|Data|Type|Description|Expected Values|
|----|----|----|----|
|Age|Int|Age of agent|10-90|
|Status|Symbol|Infection status|S, I, R|
|Severity|Float|Current infection severity|0-100|
|Mask|Bool|Uses mask?|Y,N|
|<img src="https://render.githubusercontent.com/render/math?math=\beta">|Float|Current infection Probability|0-1|
|<img src="https://render.githubusercontent.com/render/math?math=r">|Float|Reinfection probability|0-1|
|Quarantined|Bool|Currently Quarantined|Y,N|
|Hospitalized|Bool|Currently Hospitalized|Y,N|
|Speed|Float|Movement speed (dependent on global quarantine level)|0-1|
|Position|Tuple(Float, 2)|Current location|(0-1, 0-1)|
|Mass|Float|Interaction prevention parameter|0, <img src="https://render.githubusercontent.com/render/math?math=\infinity">|

## Model Steps
In each step of the model, the following global actions will be performed:
1. Check for interactions between agents
2. Transmit the virus if necessary
3. Update the quarantine level

In each step of the model, the following agent-level actions will be performed:
- Move the agent to a new position
- Update the agent's values
  - If infected, then do the following:
    - Increment infection time counter
    - Update infection severity (how sick is the patient right now?)
    - Check to see if infection is detected
    - Check to see if agent has died
    - Check to see if agent has recovered
  - If recovered, then do the following:
    - update reinfection probability
  - If not infected, then do nothing
  - Update mass value (<img src="https://render.githubusercontent.com/render/math?math=\infinity"> if quarantined or hospitalized, 1 otherwise)
  - Update speed value (slow release of quarantine)

## Transmission probability
The viral transmission probability curve will be a piecewise-defined function with 5 parameters.  These are:
1. <img src="https://render.githubusercontent.com/render/math?math=t_{\beta_{max}}">: The time when the maximum transmission probability is reached (2-6 days, based on known data)
2. <img src="https://render.githubusercontent.com/render/math?math=t_{\beta_{med}}">: The time when the median transmission probability is reached. (halfway between 2 days and <img src="https://render.githubusercontent.com/render/math?math=t_{\beta_{max}}">)
3. <img src="https://render.githubusercontent.com/render/math?math=K_\beta">: The slope of transmission probability increase (arbitrarily 5-10)
4. <img src="https://render.githubusercontent.com/render/math?math=\beta_{max}">: The maximum transmission probability (around 5%)
5. <img src="https://render.githubusercontent.com/render/math?math=\eta">: The rate of transmission probability decrease. (set so that it is near 0 by 11 days)

The probability curve will sketch out a logistic function until time
<img src="https://render.githubusercontent.com/render/math?math=t_{\beta_{max}}">, after which it will undergo a linear decrease with slope <img src="https://render.githubusercontent.com/render/math?math=\eta">.  It is mathematically defined as

<img src="http://www.sciweavers.org/tex2img.php?eq=%20%5Cbeta_t%20%3D%5Cbegin%7Bcases%7D%5Cbeta_%7Bt-1%7D%2B%5CDelta%20t%20%5Cfrac%7Be%5E%7B-%28t-t_%7B%5Cbeta_%7Bmax%7D%7D%29%2FK_%5Cbeta%7D%7D%7BK_%5Cbeta%281%2Be%5E%7B-%28t-t_%7B%5Cbeta_%7Bmax%7D%7D%29%2FK_%5Cbeta%7D%29%5E2%7D%20%26%20t%20%3C%20t_%7B%5Cbeta_%7Bmax%7D%7D%5C%5C%0A%5Cbeta_%7Bt-1%7D%20-%20%5Ceta%5CDelta%20t%20%26%20t%20%5Cge%20t_%7B%5Cbeta_%7Bmax%7D%7D%5Cend%7Bcases%7D%20&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt=" \beta_t =\begin{cases}\beta_{t-1}+\Delta t \frac{e^{-(t-t_{\beta_{max}})/K_\beta}}{K_\beta(1+e^{-(t-t_{\beta_{max}})/K_\beta})^2} & t < t_{\beta_{max}}\\\beta_{t-1} - \eta\Delta t & t \ge t_{\beta_{max}}\end{cases} " width="361" height="81" />

The increase and decrease of the transmission probability is not affected by any other status.

## Disease Severity
The disease severity curve will be a piecewise-defined function with 5 parameters.  These are:
1. <img src="https://render.githubusercontent.com/render/math?math=t_{\S_{max}}">: The time of maximum disease severity (<img src="https://render.githubusercontent.com/render/math?math=t_{S_0}"> + 1-14 days)
2. <img src="https://render.githubusercontent.com/render/math?math=t_{S_0}">: Time of the onset of noticable symptoms (2-14 days)
3. <img src="https://render.githubusercontent.com/render/math?math=K_S">: How quickly symptoms manifest and become worse (arbitrarily 0.5-3)
4. <img src="https://render.githubusercontent.com/render/math?math=S_{max}">: The maximum theoretical disease severity for this agent (based on age alone.)
5. <img src="https://render.githubusercontent.com/render/math?math=\gamma">: The rate of disease recovery.

The disease severity curve will sketch out a logistic function until time
<img src="https://render.githubusercontent.com/render/math?math=t_{S_0}">, after which it will undergo a linear decrease with slope <img src="https://render.githubusercontent.com/render/math?math=\gamma">.  It is mathematically defined as

<img src="http://www.sciweavers.org/tex2img.php?eq=S_t%20%3D%5Cbegin%7Bcases%7DS_%7Bt-1%7D%2B%5CDelta%20t%20%5Cfrac%7Be%5E%7B-%28t-t_%7BS_0%7D%29%2FK_S%7D%7D%7BK_S%281%2Be%5E%7B-%28t-t_%7BS_%7Bmax%7D%7D%29%2FK_S%7D%29%5E2%7D%20%26%20t%20%3C%20t_%7BS_%7Bmax%7D%7D%5C%5C%0AS_%7Bt-1%7D%20-%20%5Cgamma%5CDelta%20t%20%26%20t%20%5Cge%20t_%7BS_%7Bmax%7D%7D%5Cend%7Bcases%7D%20&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0" align="center" border="0" alt="S_t =\begin{cases}S_{t-1}+\Delta t \frac{e^{-(t-t_{S_0})/K_S}}{K_S(1+e^{-(t-t_{S_{max}})/K_S})^2} & t < t_{S_{max}}\\S_{t-1} - \gamma\Delta t & t \ge t_{S_{max}}\end{cases} " width="360" height="81" />

If the agent is hospitalized, then the increase in severity will be cut by 1/4.

## Reinfection
The probability of reinfection will slowly increase up to the original maximum transmission probability, where it will be half of the original probability after 30 days.

## Wearing masks
Each mask involved in the potential transmission will decrease the transmission probability by 1/3.

## Detection, hospitalization, and death
An individual will be detected as being infected if <img src="https://render.githubusercontent.com/render/math?math=S_t \ge S_{detect}"> or
<img src="https://render.githubusercontent.com/render/math?math=S_t \ge rand()">.

An individual will be hospitalized if <img src="https://render.githubusercontent.com/render/math?math=S_t \ge S_{crit}"> AND <img src="https://render.githubusercontent.com/render/math?math=rand() \le S_{t}/4">

An individual will be considered to have died from COVID if <img src="https://render.githubusercontent.com/render/math?math=S_t \ge S_{crit}"> for 5 days.




