# Schemas

Schemas define the structure, data types, and additional metadata for datasets used within a defined project space.
They ensure consistency and provide clear guidelines on how data is formatted and interpreted.

## Information

All schema items are designed to describe properties related to a single structure. Each item in the schema contains the following fields.

### Description

This multi-line text field provides a brief yet comprehensive description of the data item.
It explains what the data represents, how it's used, and any specific details relevant to understanding the significance of the data within the context of the dataset.

### Shape

The `shape` field describes the dimensions of the data with the assumption that we have multiple structures or time points.
This field helps users understand how data is organized, whether it's a single value, a vector, a matrix, or a higher-dimensional array.

We use `n_structures` and `n_atoms` as placeholders to account for an arbitrary number of structures or atoms in specific systems of the same number of atoms.
These values are a list of the aforementioned placeholders or integers for the shape of that dimension.
For example, the shape of atomic Cartesian coordinates would be `[n_structures, n_atoms, 3]` and the atomic numbers of all atoms in the system would be `[n_atoms]`.

### Dtype

The `dtype` (data type) field specifies the type of data stored in the dataset and is detailed [below](#data-types).
The data type is crucial for understanding what kind of operations can be performed on the data and how it should be handled during processing.

### Units

The `units` field indicates the measurement units for the data item. This is essential for interpreting the data correctly, especially in scientific calculations where dimensions must be consistent. Examples of units include `meters` for length, `seconds` for time, `Kelvin` for temperature, etc.
If the data does not have units (e.g., it's a count or a dimensionless quantity), this should be set to `null`. This specification helps prevent errors in data interpretation by clearly stating when units are not applicable.

In summary, these schema definitions serve as a blueprint for the data, providing all necessary details to use, interpret, and process the data correctly. Each aspect of the schema plays a critical role in data management and utilization.

## Data types

### Integers

| Type | Alias | Description |
| ---- | ----- | ----------- |
| `int8` | `i1` | An 8-bit signed integer whose values exist on the interval `[-128, +127]`. |
| `int16` | `i2` | A 16-bit signed integer whose values exist on the interval `[−32,767, +32,767]`. |
| `int32` | `i4` | A 32-bit signed integer whose values exist on the interval `[−2,147,483,647, +2,147,483,647]`. |
| `int64` | `i8` | A 64-bit signed integer whose values exist on the interval `[−9,223,372,036,854,775,807, +9,223,372,036,854,775,807]`. |

| Type | Alias | Description |
| ---- | ----- | ----------- |
| `uint8` | `u1` | An 8-bit unsigned integer whose values exist on the interval `[0, +255]`. |
| `uint16` | `u2` | A 16-bit unsigned integer whose values exist on the interval `[0, +65,535]`. |
| `uint32` | `u4` | A 32-bit unsigned integer whose values exist on the interval `[0, +4,294,967,295]`. |
| `uint64` | `u8` | A 64-bit unsigned integer whose values exist on the interval `[0, +18,446,744,073,709,551,615]`. |

### Floats

| Type | Alias | Description | Min > 0 | Max > 0 |
| ---- | ----- | ----------- | ------- | ------- |
| `float32` | `f4` | IEEE 754 single-precision binary floating-point number. | 1.18 ⨉ 10<sup>-38</sup> | 3.40 ⨉ 10<sup>38</sup> |
| `float64` | `f8` | IEEE 754 double-precision binary floating-point number. | 2.23 ⨉ 10<sup>-308</sup> | 1.80 ⨉ 10<sup>308</sup> |

#### Precision

The figure and table below shows the absolute precision for both formats over a range of values given by

$$
\text{precision} = \frac{\text{value}}{2^{n}};
$$

where $n$ is the number of mantissa bits which is 23 for `float32` and 52 for `float64`.
This can be used to select an appropriate format given the expected value of a number and the required precision.

| Value | `float32` precision | `float64` precision |
| ----- | ------------------- | --------------- |
|  **10<sup>-12</sup>** | 10<sup>-19</sup> | 10<sup>-28</sup> |
|  **10<sup>-9</sup>** | 10<sup>-16</sup> | 10<sup>-25</sup> |
|  **10<sup>-6</sup>** | 10<sup>-13</sup> | 10<sup>-22</sup> |
|  **10<sup>-3</sup>** | 10<sup>-10</sup> | 10<sup>-19</sup> |
|  **10<sup>0</sup>** | 10<sup>-7</sup> | 10<sup>-16</sup> |
|  **10<sup>3</sup>** | 10<sup>-4</sup> | 10<sup>-13</sup> |
|  **10<sup>6</sup>** | 10<sup>-1</sup> | 10<sup>-10</sup> |
|  **10<sup>9</sup>** | 10<sup>2</sup> | 10<sup>-7</sup> |
|  **10<sup>12</sup>** | 10<sup>5</sup> | 10<sup>-4</sup> |

<div id="container-float-precision"></div>
<script type="module">
// Set up the dimensions for the SVG container
const width = 600;
const height = 600;
const margin = { top: 20, right: 20, bottom: 50, left: 60 };
// Calculate the actual width and height available for the chart
const innerWidth = width - margin.left - margin.right;
const innerHeight = height - margin.top - margin.bottom;
// Create the SVG container
const svg = d3.select('#container-float-precision').append('svg')
  .attr('width', width)
  .attr('height', height);
// Create a group element to contain the chart and apply margins
const chartGroup = svg.append('g')
  .attr('transform', `translate(${margin.left},${margin.top})`);
// Set up scales for x and y axes
const xScale = d3.scaleLog()
  .domain([10**-12, 10**12]) // Assumes the x-values are the same for both lines
  .range([0, innerWidth])
  .base(10);
const yScale = d3.scaleLog()
  .domain([10**-29, 10**6])
  .range([innerHeight, 0])
  .base(10);
// Create line generator function
const line = d3.line()
  .x(d => xScale(d.x))
  .y(d => yScale(d.y));
// Draw x-axis
chartGroup.append('g')
  .attr('transform', `translate(0,${innerHeight})`)
  .call(d3.axisBottom(xScale).tickFormat(d3.format(".0e")));
chartGroup.append("text")
    .attr("text-anchor", "end")
    .attr("y", 560)
    .attr("x", 350)
    .attr("dy", ".75em")
    .text("Floating point value");
// Draw y-axis
chartGroup.append('g')
  .call(d3.axisLeft(yScale).tickFormat(d3.format(".0e")));
chartGroup.append("text")
    .attr("text-anchor", "end")
    .attr("y", -60)
    .attr("x", -180)
    .attr("dy", ".75em")
    .attr("transform", "rotate(-90)")
    .text("Floating point precision");
// Draw lines
// float32
chartGroup.append('path')
  .datum([{ x: 10**-12, y: 1.1920928955078125*10**-19 }, { x: 10**12, y: 119209.28955078125 }])
  .attr('fill', 'none')
  .attr('stroke', `#e41a1c`)
  .attr('stroke-width', 2)
  .attr('d', line);
chartGroup.append("text")
  .attr("transform", `translate(60,280)`)
  .attr("text-anchor", "start")
  .style("fill", "#e41a1c")
  .text("float32");
// float64
chartGroup.append('path')
  .datum([{ x: 10**-12, y: 2.220446049250313*10**-28 }, { x: 10**12, y: 0.0002220446049250313 }])
  .attr('fill', 'none')
  .attr('stroke', `#377eb8`)
  .attr('stroke-width', 2)
  .attr('d', line);
chartGroup.append("text")
  .attr("transform", `translate(60,420)`)
  .attr("text-anchor", "start")
  .style("fill", "#377eb8")
  .text("float64");
</script>

### Text

We generally use Unicode with type `utf8`.

### Boolean

Type `bool` is used for `True` (`1`) and `False` (`0`).
