# Healthcare Innovation Journey

My name is **Eduardo Corazin**, and this project was built as part of a technical challenge for **Mastera**.  
It is an interactive dashboard designed to guide users (researchers, doctors, pharma companies) in exploring opportunities for medical innovation through scientific articles and patents.

## Project Overview

An interactive **Streamlit app** built in Python that integrates public data sources (**PubMed API + semi-realistic patents**) and stores user interactions in a relational database (**PostgreSQL on Supabase**).

### Note on Patents
The patent dataset in this project is **semi-realistic and simulated**.  
This decision was made because APIs such as **PatentsView** require key approval, which usually takes 1–3 business days.  
Given the short time frame of the challenge, I implemented a generator that creates realistic-looking patents based on article keywords.  
This module can be easily replaced by a direct PatentsView API integration once access keys are available.

---

## Features

- Interactive search: enter a healthcare topic (e.g., "diabetes treatment") and choose focus (Treatment, Diagnosis, Prevention).  
- Scientific articles: collected via [PubMed API](https://www.ncbi.nlm.nih.gov/pubmed/), automatically translated to Portuguese.  
- Patent suggestions: semi-realistic simulated patents derived from article keywords.
- Keyword Extraction: identify key concepts from article titles using KeyBERT.
- Gap Analysis: highlights opportunities where academic research is high but patent coverage is low.
- Dashboard: visual metrics, bar charts (Plotly) and word clouds for keyword frequency.  
- Detailed tables: interactive data grids with sorting and filtering (AgGrid).  
- Report Generation: export insights and recommendations as TXT files.
- Database logging: all user interactions (queries + counts of articles/patents) are anonymized and stored in PostgreSQL (Supabase free-tier).  

---

## Tech Stack

- **Frontend & UI:** Streamlit  
- **Data Processing:** pandas, requests, biopython (Entrez), googletrans
- **NLP & Keyword Extraction:** keybert, sentence-transformers, torch
- **Visualizations:** plotly, wordcloud, streamlit-aggrid  
- **Database:** PostgreSQL (Supabase) + SQLAlchemy  
- **Environment management:** python-dotenv  

---

## Database

The project uses **Supabase (PostgreSQL)** as the database provider.  
A table `user_interactions` must exist with the following schema:

```sql
CREATE TABLE user_interactions (
    id SERIAL PRIMARY KEY,
    user_query TEXT,
    anonymized_id VARCHAR(64),
    results_count INT,
    data_sources JSONB,
    created_at TIMESTAMP DEFAULT NOW()
);
````
Queries are anonymized with SHA-256 hashing.
Only aggregated counts are stored (compliant with LGPD/GDPR).

Setup & Installation
Clone the repo and install dependencies:

```git
git clone https://github.com/eduardo-corazin/jornada-inovacao-saude.git
cd jornada-inovacao-saude
pip install -r requirements.txt
````
Create a .env file (never commit this) with:

```git
SUPABASE_USER=your_user
SUPABASE_PASS=your_password
SUPABASE_HOST=your_host.supabase.co
SUPABASE_PORT=5553
SUPABASE_DB=your_db
NCBI_EMAIL=your_email@example.com
USER_AGENT=healthcare-innovation
````
Run locally:

```vs
streamlit run src/app/app.py
````
# Deployment

The app is deployed on Streamlit Cloud

# Challenge Context

This project was developed for the Mastera Internship Challenge.
It simulates a real-world tool used by pharmaceutical companies, researchers, and hospitals to explore innovation opportunities in healthcare.

Interactive journeys: guided exploration of health topics

External data sources: PubMed + semi-realistic patents (replaceable with PatentsView API)

Database logging: user queries stored securely

Visualization & insights: dashboards and downloadable reports

# Challenge Context

**Interactive journeys: guided exploration of health topics**

**External data sources: PubMed + semi-realistic patents (replaceable with PatentsView API)**

**Database logging: user queries stored securely**

**Visualization & insights: dashboards, heatmaps, word clouds, and downloadable reports**

**Gap analysis: detect innovation opportunities where academic research exceeds patent coverage**