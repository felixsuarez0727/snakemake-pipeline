# Use an official Python runtime as a parent image
FROM python:3.9

# Set the working directory in the container
WORKDIR /usr/src/app

# Install BWA and Samtools
RUN apt-get update && \
    apt-get install -y bwa samtools && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install the pulp library for linear programming
RUN pip install --no-cache-dir pulp==2.4


# Copy the current directory contents into the container at /usr/src/app
COPY . /usr/src/app

# Define environment variable for Snakemake
ENV LC_ALL=C.UTF-8

# Install Snakemake
RUN pip install --no-cache-dir snakemake

# Run Snakemake pipeline by default when the container starts
CMD ["snakemake", "--cores", "1", "-s", "Snakefile"]
