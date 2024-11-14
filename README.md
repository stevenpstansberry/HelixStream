# HelixStream: A Bioinformatics Pipeline for Analyzing Viral Evolution

## Overview

Welcome to **HelixStream**! HelixStream is a bioinformatics pipeline designed to analyze genetic mutations and evolutionary patterns in viral genomes, with a particular focus on SARS-CoV-2 variants such as Omicron. Built to be flexible and extensible, HelixStream can easily adapt to other bioinformatics use cases. By leveraging modern technologies across cloud infrastructure, data processing, and interactive visualization, HelixStream offers a comprehensive solution for tracking, analyzing, and visualizing genetic sequence data. Whether you’re a researcher or a data enthusiast, HelixStream empowers you to observe viral evolution trends, track mutation rates, and compare genetic distances across different variants in real-time.

## Table of Contents

- [Features](#features)
- [Tech Stack](#tech-stack)
- [Setup and Installation](#setup-and-installation)
- [Future Improvements](#future-improvements)

## Features

### 1. Variant Mutation Analysis 🔍

HelixStream tracks the mutations in viral genomes, focusing on specific SARS-CoV-2 variants. The pipeline compares sequences over time, enabling users to observe evolutionary trends, mutation rates, and key genetic differences across variants like Omicron, Delta, and Alpha. 

### 2. Real-Time Sequence Ingestion and Analysis 🚀

Data is fetched directly from the **NCBI** database via Linux/Unix scripts running on an **AWS EC2** instance. These scripts automate the ingestion of new viral sequences into an **S3 bucket** for scalable storage, providing a seamless flow of data from source to analysis.

### 3. Sequence Alignment and Mutation Metrics Calculation 🔬

Using **Python** and **Biopython**, HelixStream globally aligns sequences through **Clustal Omega** and computes essential metrics such as **Jukes-Cantor Distance** for evolutionary divergence, mutation density per kilobase, and percentage sequence identity. The data is cleaned and preprocessed, making it ready for downstream analyses.

### 4. API Backend with Spring Boot 🌐

HelixStream’s API backend, built with **Spring Boot**, acts as a bridge to execute Linux commands on the EC2 instance, trigger data ingestion, and handle processing requests. This backend service provides a robust API layer that can be used to interact programmatically with the HelixStream pipeline.

### 5. Interactive Data Visualization Dashboard 📊

HelixStream includes a **React** frontend for real-time visualization of key metrics, enabling users to track metrics like Percent Identity, Mutation Density, and Jukes-Cantor Distance over time. With multi-variant comparison capabilities, users can switch between different COVID-19 variants to compare mutation profiles, making it a powerful tool for observing viral evolution.

<p align="center">
  <img src="./images/helixstream-dashboard.png" alt="HelixStream Dashboard" width="70%">
</p>

## Tech Stack

### Frontend

- **React**: Provides an interactive user interface for selecting variants, viewing mutation metrics, and comparing different datasets.
- **Next.js**: Seamless integration of backend logic for dynamic bioinformatics data visualization and variant comparison.

### Backend

- **Spring Boot**: Serves as the API layer, managing commands and interactions with the EC2 instance for data ingestion and processing tasks.
- **Linux/Unix Scripting**: Used on the EC2 instance to fetch data from NCBI, upload sequences to S3, and perform other processing tasks.

### Data Processing

- **Python**: Handles data cleaning, sequence alignment, and mutation metric calculation.
- **Biopython**: Provides tools for sequence alignment and analysis, specifically with Clustal Omega, to calculate evolutionary distances.

### Infrastructure

- **AWS S3**: Stores ingested and processed data in a scalable, durable way.
- **AWS EC2**: Runs bash scripts to fetch data, execute processing jobs, and upload results to S3.


## Setup and Installation

### Prerequisites

- **AWS Account** with IAM privileges for EC2, S3
- **React** installed on your machine.
- **AWS CLI** to configure credentials.

### Clone the Repository

```bash
git clone https://github.com/stevenpstansberry/HelixStream.git
cd HelixStream/Frontend
npm install
```

Run the Application
Start the React development server:

```bash
npm start
This will run the app locally at http://localhost:3000 by default.
```

Run Backend Services
Deploy the Spring Boot backend to handle API requests and EC2 scripts for data ingestion and processing.

## Future Improvements

- **Expanded Variant Support**: Add support for tracking additional viral or bacterial pathogens.

- **Automated Data Refresh**: Implement scheduled data ingestion from NCBI, allowing for real-time updates without manual intervention.

- **Advanced Analytics**: Integrate machine learning algorithms for predictive modeling of mutation trends and variant evolution.

- **Public Dashboard Access**: Enhance Grafana settings or use a self-hosted solution for sharing dashboards publicly.

---

**HelixStream** represents a powerful, extensible pipeline for bioinformatics research and analytics, making it adaptable for any project requiring sequence alignment, evolutionary distance metrics, or multi-variant comparison. Its design allows for future expansion, making it a robust tool for tracking viral evolution or applying to other bioinformatics domains.

