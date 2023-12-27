# Schemas

TODO:

## Information

All schema items contain the following information to describe a single structure.

### description

TODO:

### ndim

TODO:

### dtype

TODO:

### units

TODO:

## Data types

### Integers

| Type | Alias | Description |
| ---- | ----- | ----------- |
| `i1` | `int8` | An 8-bit signed integer whose values exist on the interval `[-128, +127]`. |
| `i2` | `int16` | A 16-bit signed integer whose values exist on the interval `[−32,767, +32,767]`. |
| `i3` | `int32` | A 32-bit signed integer whose values exist on the interval `[−2,147,483,647, +2,147,483,647]`. |
| `i4` | `int64` | A 64-bit signed integer whose values exist on the interval `[−9,223,372,036,854,775,807, +9,223,372,036,854,775,807]`. |

| Type | Alias | Description |
| ---- | ----- | ----------- |
| `u1` | `uint8` | An 8-bit unsigned integer whose values exist on the interval `[0, +255]`. |
| `u2` | `uint16` | A 16-bit unsigned integer whose values exist on the interval `[0, +65,535]`. |
| `u4` | `uint32` | A 32-bit unsigned integer whose values exist on the interval `[0, +4,294,967,295]`. |
| `u8` | `uint64` | A 64-bit unsigned integer whose values exist on the interval `[0, +18,446,744,073,709,551,615]`. |

### Floats

| Type | Alias | Description | Min > 0 | Max > 0 |
| ---- | ----- | ----------- | ------- | ------- |
| `f4` | `float32` | IEEE 754 single-precision binary floating-point number. | 1.18 ⨉ 10<sup>-38</sup> | 3.40 ⨉ 10<sup>38</sup> |
| `f8` | `float64` | IEEE 754 double-precision binary floating-point number. | 2.23 ⨉ 10<sup>-308</sup> | 1.80 ⨉ 10<sup>308</sup> |

#### Precision

The figure and table below shows the absolute precision for both formats over a range of values given by

$$
\text{precision} = \frac{\text{value}}{2^{n}};
$$

where $n$ is the number of mantissa bits which is 23 for `f4` and 52 for `f8`.
This can be used to select an appropriate format given the expected value of a number and the required precision.

| Value | `f4` precision | `f8` precision |
| ----- | ------------- | --------------- |
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
  .attr("transform", `translate(100,280)`)
  .attr("text-anchor", "start")
  .style("fill", "#e41a1c")
  .text("f4");
// float64
chartGroup.append('path')
  .datum([{ x: 10**-12, y: 2.220446049250313*10**-28 }, { x: 10**12, y: 0.0002220446049250313 }])
  .attr('fill', 'none')
  .attr('stroke', `#377eb8`)
  .attr('stroke-width', 2)
  .attr('d', line);
chartGroup.append("text")
  .attr("transform", `translate(100,420)`)
  .attr("text-anchor", "start")
  .style("fill", "#377eb8")
  .text("f8");
</script>

### Text

We generally use Unicode with type `U` with a fixed length.
For example, Unicode text with a maximum of 30 characters would be `U30`.

### Boolean

Type `b` is used for `True` (`1`) and `False` (`0`).
