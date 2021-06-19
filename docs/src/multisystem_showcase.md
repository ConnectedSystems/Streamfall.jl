# A simple showcase of multi-system considerations

This page is a draft.

Here we showcase a two-node network representing a river and a dam downstream.

The Lower Campaspe catchment - a small semi-arid basin in North-Central Victoria, Australia - is used for the example here.

- Figure of catchment

As a graph, the network looks like this:

- Figure of two-node network

The dam is the primary water store for farmers in the area but is also used for recreational activities (camping, boating, fishing, etc) by local enthusiasts and vacationers. The Campaspe river is also home to a culturally and ecologically significant population of fish. A certain level of flow must be ensured at key times during the year to support and maintain their population levels.

In this hypothetical study, local stakeholders would like to have an idea of the range of possible
dam levels under a range of environmental watering policies, farmer water use, and how this may impact the level of enjoyment by vacationers.


The possible environmental watering strategies are defined as:

- Implicit watering:
    No purposeful releases for environmental demands. Assume natural inflows and agricultural water orders provide sufficient water flow for ecological purposes.

- Explicit watering:
    Assume agricultural water orders partially fulfill environmental needs. Water is released as needed to meet any deficit.

- Prioritized watering:
    Water for environmental purposes are prioritised and are managed separately from agricultural demands.


For the purpose of this example, the farm water requirements are given as a volume of daily water releases throughout a growing season; the period of time over which a crop can grow. This figure may be provided by another model in practice. The growing season is assumed to be between *X* and *Y*, with the daily water requirements over that period being between *X* and *Y*.

An index value is used to provide indications of the suitability of dam levels for recreational purposes.

- explain how recreational index works

Another indicator model is used to show how often environmental needs are met.

- explain how the environmental indicator model works

An overview of the system under investigation can then be conceptualized like:

- Conceptual figure of the system

Where water flows into the dam, and water is released to fulfill water needs of the users downstream. Note that "water users" as defined here includes the environment itself.

First, we define a two-node Streamfall Network which represents the river and dam:

```yaml
```

We can then generate a number of scenarios representing a mix of the management strategies and water demands, as listed above.

```julia
# Code to generate scenarios
```

```julia
# code showing how to run the model(s)
```

Analysis and wrap up...