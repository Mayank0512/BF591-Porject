#!/bin/bash

# Create package directory
package_dir="home/mayank/BU/NonLab/BF591/project/helper_packages"
mkdir "$package_dir"
cd "$package_dir"

# Create R directory
mkdir R

# Create R functions file
touch R/functions.R

# Create man directory
mkdir man

# Create man files
touch man/meta_info_from_labels.Rd
touch man/sample_replicate.Rd
touch man/timepoint_from_sample.Rd

# Create DESCRIPTION file
echo "Package: meta_data_from_sample
Type: Package
Title: Sample Information Package
Version: 0.1.0
Author:Mayank Ghogale
Maintainer: maubu@bu.edu
Description: A package to extract meta-data from sample names." > DESCRIPTION

# Create NAMESPACE file
echo "export(timepoint_from_sample)
export(sample_replicate)
export(meta_info_from_labels)" > NAMESPACE

echo "Package structure created successfully!"

