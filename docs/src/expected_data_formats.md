# Input data format

By convention, Streamfall expects data to be formatted in the following manner:

- All time series must have a "Date" column in `YYYY-mm-dd` format.
- By convention, column names follow the pattern: `[node_name]_[phenomenon]_[other metadata]`
- Unit of measure itself is optionally included in square brackets (`[]`)
- In-file comments are indicated with a hash (`#`)

## Streamflow data

Data for a node representing a subcatchment or gauge should, at a minimum, include columns
for date and streamflow (`Q`) in a consistent unit (in this case, megaliters/day).

Streamflow columns should begin with its name and include an identifiable suffix.
Here, `_Q` is used but this can be user-defined when preparing data for modelling.

```csv
# A comment
Date,406219_Q_[ML]
1988-07-23,0.5184
1988-07-24,13.6512
1988-07-25,407.6352
1988-07-26,2488.4928
1988-07-27,3443.6448
... snip ...
```

## Dam releases

Note that releases is denoted with `_releases`.

```csv
# A comment
Date, 406000_releases_[ML], 406000_extractions_[ML]
1981-01-01,0,0
1981-01-02,359,0
1981-01-03,359,0
1981-01-04,359,0
1981-01-05,460,0
... snip ...
```

## Climate data

Streamfall currently expects data for all nodes to be stored in one large DataFrame.


### Example climate data structure

```csv
Date, 406214_P, 406214_ET, 406219_P, 406219_ET
1981-01-01, 0.0, 4.8, 0.0, 4.9
1981-01-02, 0.1, 0.5, 0.1, 3.3
1981-01-03, 10.5, 5.3, 7.2, 2.3
1981-01-04, 9.89, 7.9, 6.1, 4.3
1981-01-05, 0.3, 4.2, 0.2, 6.4
... snip ...
```

where `P` indicates rainfall, and `ET` denotes evapotranspiration.

See [Climate](@ref) for more detail on climate data functions.
