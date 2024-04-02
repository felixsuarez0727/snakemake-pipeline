# Use an official Python runtime as a parent image
FROM python:3.9

# Set the working directory in the container
WORKDIR /usr/src/app

# Install BWA version 0.7.17
RUN apt-get update && \
    apt-get install -y bwa=0.7.17* && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Download and install Samtools version 1.19.2
RUN wget https://github.com/samtools/samtools/releases/download/1.19.2/samtools-1.19.2.tar.bz2 && \
    tar -xvf samtools-1.19.2.tar.bz2 && \
    cd samtools-1.19.2 && \
    ./configure && \
    make && \
    make install && \
    cd .. && \
    rm -rf samtools-1.19.2 samtools-1.19.2.tar.bz2

# Install the pulp library for linear programming version 2.4
RUN pip install --no-cache-dir pulp==2.4

# Copy the current directory contents into the container at /usr/src/app
COPY . /usr/src/app

# Define environment variable for Snakemake
ENV LC_ALL=C.UTF-8

# Install Snakemake version 6.0.0
RUN pip install --no-cache-dir snakemake==6.0.0

# Run Snakemake pipeline by default when the container starts
CMD ["snakemake", "--cores", "1", "-s", "Snakefile"]
