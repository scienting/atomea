# Digesters

In the field of computational chemistry and molecular simulations, data is often stored in a variety of file formats, each with its own structure and conventions.
This heterogeneity can lead to challenges in data analysis, as researchers must write custom code to handle each format.
To address this issue, we introduce the concept of "Digesters"&mdash;a modular and abstract way to convert data from various formats into a standardized, compressed format that can be easily analyzed using a common set of tools.

## What are Digesters?

Digesters are modular components designed to convert data from a specific input format into a standardized output format.
Each digester is tailored to handle a particular file format commonly used in computational chemistry or molecular simulations, such as XYZ, PDB, Amber simulations, etc.
The digester reads the input file(s), extracts the relevant information, and transforms it into a file format requested by the user.

By encapsulating the conversion process within a digester, we abstract away the complexities of dealing with different file formats, allowing users to focus on the analysis of their data rather than the intricacies of parsing each format.

## Benefits of using Digesters

### Standardization

Digesters ensure that all data, regardless of its original format, is converted into a consistent, standardized format. This standardization simplifies the development of analysis tools, as they only need to handle a single data format.

### Modularity

Each Digester is a self-contained module responsible for converting a specific file format. This modularity makes it easy to add support for new formats without modifying the core analysis code. Users can create their own Digesters for custom file formats, extending the functionality of the package.

### Compression

Digesters not only convert the data but also compress it to reduce storage requirements. This compression is particularly beneficial when dealing with large datasets, such as molecular simulation trajectories, which can occupy significant disk space.

### Ease of use

By abstracting away the complexities of file format conversion, Digesters provide a simple and intuitive interface for users. Researchers can focus on their scientific questions rather than spending time on data preprocessing.

## How to use Digesters

Using Digesters is straightforward.
First, identify the format of your input data and select the corresponding Digester.
Then, create an instance of the Digester, specifying the input file path and any additional parameters required by the specific Digester.

Here's a general example of how to use a digester:

TODO: Add example

Each Digester may have its own set of parameters and options, so be sure to consult the documentation for the specific Digester you are using.

## Creating custom Digesters

If you have data in a custom format that is not yet supported by the built-in digesters, you can create your own to handle the conversion process.
To create a custom Digester, follow these steps:

TODO:
