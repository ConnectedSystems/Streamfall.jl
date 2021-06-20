# Input data format

By convention, Streamfall expects data to be formatted in the following manner:

- All time series must have a "Date" column in `YYYY-mm-dd` format.
- By convention, column names follow the pattern: `[node_name]_[phenomenon]_[other metadata]`
- Unit of measure itself is optionally included in square brackets (`[]`)
- In-file comments are indicated with a hash (`#`)


## Example

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
Date, 406214_rain, 406214_evap, 406219_rain, 406219_evap
1981-01-01, 0.0, 4.8, 0.0, 4.9
1981-01-02, 0.1, 0.5, 0.1, 3.3
1981-01-03, 10.5, 5.3, 7.2, 2.3
1981-01-04, 9.89, 7.9, 6.1, 4.3
1981-01-05, 0.3, 4.2, 0.2, 6.4
... snip ...
```

### Climate data functions

```@autodocs
Modules = [Streamfall]
Order   = [:function, :type]
Pages   = ["Climate.jl"]
```