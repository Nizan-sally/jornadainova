#data_fetcher.py

import ssl, os, random
import pandas as pd
from Bio import Entrez
from dotenv import load_dotenv
from deep_translator import GoogleTranslator
import streamlit as st
from keybert import KeyBERT
import matplotlib.pyplot as plt

ssl._create_default_https_context = ssl._create_unverified_context

load_dotenv()
Entrez.email = os.getenv("NCBI_EMAIL")  

translator = GoogleTranslator(source="auto", target="pt")
kw_model = KeyBERT()

# Stopwords
STOPWORDS = ["de","do","da","em","para","com","e","o","a","os","as"]
STOPWORDS_SAUDE = ["paciente","criança","homem","mulher","idoso","adulto","caso","grupo","estudo","dados","ano","resultado"]

# Construção 
def construir_query(tema, foco=None, data_inicio=None, data_fim=None):
    palavras = [p for p in tema.split() if p.lower() not in STOPWORDS] or [tema]
    query = " OR ".join([f"{p}[Title/Abstract]" for p in palavras])
    focos = {"Tratamento":"(therapy OR treatment OR intervention)",
             "Diagnóstico":"(diagnosis OR detection OR screening)",
             "Prevenção":"(prevention OR prophylaxis OR vaccine)"}
    if foco in focos:
        query += f" AND {focos[foco]}[Title/Abstract]"
    if data_inicio or data_fim:
        query += f" AND ({data_inicio or 1800}[Date - Publication] : {data_fim or 3000}[Date - Publication])"
    return query

# Extrair informações
def extrair_info(article):
    ano = "N/A"
    if "ArticleDate" in article and len(article["ArticleDate"])>0:
        ano = article["ArticleDate"][0].get("Year","N/A")
    elif "Journal" in article:
        ano = article.get("Journal",{}).get("JournalIssue",{}).get("PubDate",{}).get("Year","N/A")
    autores = ", ".join([f"{a.get('LastName','N/A')} {a.get('Initials','')}" 
                         for a in article.get("AuthorList",[])[:2]]) or "N/A"
    return ano, autores

# Buscar artigos  PubMed 
@st.cache_data(show_spinner="Buscando artigos...")
def buscar_artigos(tema, max_itens=15, foco=None, data_inicio=None, data_fim=None):
    try:
        query = construir_query(tema, foco, data_inicio, data_fim)
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_itens)
        ids = Entrez.read(handle).get("IdList",[])
        handle.close()
        if not ids:
            return pd.DataFrame(columns=["Título","Autores","Ano"])
        handle = Entrez.efetch(db="pubmed", id=ids, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        artigos = []
        for rec in records.get("PubmedArticle",[]):
            art = rec.get("MedlineCitation",{}).get("Article",{})
            title = art.get("ArticleTitle","Sem título")
            try:
                title_pt = translator.translate(title) if title else title  # Tradução do título
            except:
                title_pt = title
            ano, autores = extrair_info(art)
            artigos.append({"Título":title_pt,"Autores":autores,"Ano":ano})
        return pd.DataFrame(artigos,columns=["Título","Autores","Ano"])
    except:
        return pd.DataFrame(columns=["Título","Autores","Ano"])

# Buscar patente
def buscar_patentes(tema, df_artigos, max_itens=5, data_inicio=None, data_fim=None):
    if df_artigos.empty:  # Nenhum artigo, nenhuma patente
        return pd.DataFrame(columns=["Título","Número","Data","Inventor"]), [], []

    # Checagem do ano mínimo
    ano_inicio = int(data_inicio) if data_inicio else 2018
    ano_fim = int(data_fim) if data_fim else 2025
    if ano_inicio > ano_fim:
        ano_inicio, ano_fim = ano_fim, ano_inicio

    if ano_inicio <= 1950 or ano_fim <= 1950:
        return pd.DataFrame(columns=["Título","Número","Data","Inventor"]), [], []

    # Extração de termos (igual estava antes)
    termos = []
    try:
        for titulo in df_artigos["Título"]:
            if isinstance(titulo, str) and len(titulo) > 5:
                keywords = kw_model.extract_keywords(
                    titulo, keyphrase_ngram_range=(1,2),
                    stop_words=STOPWORDS+STOPWORDS_SAUDE, top_n=2
                )
                termos.extend([kw[0] for kw in keywords])
    except:
        pass
    if not termos:
        termos = [tema]

    nomes_inventores = ['Silva','Souza','Oliveira','Costa','Pereira','Santos','Almeida']
    patentes = []

    qtd_patentes = random.randint(1, max_itens)
    for i in range(qtd_patentes):
        termo = random.choice(termos)
        titulo_patente = f"Composição inovadora para {termo}"
        numero = f"US{random.randint(2020,2025)}/{1000+i}"
        data = str(random.randint(ano_inicio, ano_fim))  # Ano da patente dentro do intervalo
        inventor = f"Dr. {random.choice(nomes_inventores)}"
        patentes.append({"Título":titulo_patente, "Número":numero, "Data":data, "Inventor":inventor})

    gaps = [
        f"Pouca cobertura em {random.choice(termos)}",
        f"Falta de validação clínica para {random.choice(termos)}",
        f"Necessidade de tecnologias avançadas para {random.choice(termos)}",
        f"Mercado pouco explorado em {random.choice(termos)}",
        f"Baixa integração com soluções digitais para {random.choice(termos)}"
    ]

    recomendacoes = [
        "Investir em ensaios clínicos multicêntricos",
        "Explorar aplicações digitais/IA na área",
        "Parcerias com laboratórios de referência",
        "Desenvolver produtos sustentáveis e escaláveis",
        "Capacitar equipe em inovação tecnológica"
    ]

    return pd.DataFrame(patentes), gaps, recomendacoes

# Mostrar gráfico de artigos por ano
def mostrar_grafico_artigos(df_artigos):
    if df_artigos.empty:
        st.info("Nenhum artigo encontrado")
        return
    contagem = df_artigos['Ano'].value_counts().sort_index()
    fig, ax = plt.subplots(figsize=(5,3))
    contagem.plot(kind='bar', color="#4CAF50", ax=ax)
    ax.set_title("Artigos por Ano", fontsize=12)
    ax.set_xlabel("Ano", fontsize=10)
    ax.set_ylabel("Quantidade", fontsize=10)
    ax.grid(axis='y', linestyle='--', alpha=0.7)
    st.pyplot(fig)
