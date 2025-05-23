% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{capirise}
\alias{capirise}
\alias{capirise_HBV}
\alias{capirise_HBVfix}
\alias{capirise_AcceptRatio}
\alias{capirise_AcceptPow}
\title{**capilarise**}
\usage{
capirise_HBV(
  GROUND_water_mm,
  SOIL_water_mm,
  SOIL_capacity_mm,
  SOIL_potentialCapirise_mm
)

capirise_HBVfix(
  GROUND_water_mm,
  SOIL_water_mm,
  SOIL_capacity_mm,
  SOIL_fieldCapacityPerc_1,
  SOIL_potentialCapirise_mm
)

capirise_AcceptRatio(
  GROUND_water_mm,
  SOIL_water_mm,
  SOIL_capacity_mm,
  SOIL_fieldCapacityPerc_1,
  param_CAPIRISE_acr_k
)

capirise_AcceptPow(
  GROUND_water_mm,
  SOIL_water_mm,
  SOIL_capacity_mm,
  SOIL_fieldCapacityPerc_1,
  param_CAPIRISE_acp_k,
  param_CAPIRISE_acp_gamma
)
}
\arguments{
\item{GROUND_water_mm}{(mm/m2/TS) water volume in `groundLy`}

\item{SOIL_water_mm}{(mm/m2) water volume in `soilLy`}

\item{SOIL_capacity_mm}{(mm/m2) average soil Capacity (maximal storage capacity)}

\item{SOIL_potentialCapirise_mm}{<0.01, 7> (mm/m2/TS) **potential** capillary rise}

\item{SOIL_fieldCapacityPerc_1}{<0, 1> the relative ratio that the water content can drainage by gravity}

\item{param_CAPIRISE_acr_k}{<0.01, 1> coefficient parameter [capirise_AcceptRatio()]}

\item{param_CAPIRISE_acp_k}{<0.01, 1> coefficient parameter for [capirise_AcceptPow()]}

\item{param_CAPIRISE_acp_gamma}{<0.01, 1> exponential parameter for [capirise_AcceptPow()]}
}
\value{
GROUND_capirise_mm (mm/m2/TS) capillary rise
}
\description{
\loadmathjax

In hydrological modeling, capillary rise refers to the process by which water is drawn upward from groundwater (table) through the soil due to the force of capillary action.
In conceptual watershed models, the capillary rise term often refers to a process that moves water from lower to higher soil water stores, 
which may also implicitly include lateral groundwater flow processes in a sloping domain.

It can be calculated by the water in the ground layer \mjseqn{W_{grnd}}, which can also be treated as part of \mjseqn{W_{grnd}}. 
There are not many methods to describe this process. Most HMs ignore this process, 
perhaps because it is not significant in most situations, or because the process of percolation can deal with this process at the same time.

So we can give the function from:

\mjsdeqn{F_{capi} = f_{capirise}(D_{grnd}, D_{soil})}


to:

\mjsdeqn{F_{capi} = f_{capirise}(W_{grnd}, W_{soil}, C_{soil}, ...)}
\mjsdeqn{F_{capi} \leq W_{grnd}}
\mjsdeqn{F_{capi} \leq C_{soil} - W_{soil}}


where
- \mjseqn{F_{capi}} is `GROUND_capirise_mm`
- \mjseqn{W_{grnd}} is `GROUND_water_mm`
- \mjseqn{W_{soil}} is `water_SOIL_mm`
- \mjseqn{C_{soil}} is `capacity_SOIL_mm`

The output density distribution from 4 methods:
}
\details{
# **_HBV** \insertCite{HBV_Lindstrom_1997}{HydroGallery}: 


\mjsdeqn{F_{capi} = M_{capi} \left( 1 - \frac{W_{soil}}{C_{soil}} \right)}
where
  - \mjseqn{M_{capi}} is `SOIL_potentialCapirise_mm`

# **_HBVfix** \insertCite{HBV_Lindstrom_1997}{HydroGallery}: 


\mjsdeqn{F_{capi} = M_{capi} \left( 1 - \frac{W_{soil}}{k_{fc}C_{soil}} \right), \quad W_{soil} < k_{fc}C_{soil}}
where
  - \mjseqn{k_{fc}} is `SOIL_fieldCapacityPerc_1`

# **_AcceptRatio**: 


\mjsdeqn{F_{capi} = k \left( W_{soil} - k_{fc}C_{soil} \right), \quad W_{soil} < k_{fc}C_{soil}}
where
  - \mjseqn{k} is `param_CAPIRISE_acr_k`
  - \mjseqn{k_{fc}} is `SOIL_fieldCapacityPerc_1`

# **_AcceptRatio**: 


\mjsdeqn{F_{capi} = k \left( W_{soil} - k_{fc}C_{soil} \right)^\gamma, \quad W_{soil} < k_{fc}C_{soil}}
where
  - \mjseqn{k} is `param_CAPIRISE_acp_k`
  - \mjseqn{\gamma} is `param_CAPIRISE_acp_gamma`
}
\references{
\insertAllCited{}
}
