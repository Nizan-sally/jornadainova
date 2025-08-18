# app.py

import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import streamlit as st
import pandas as pd
import plotly.express as px
from wordcloud import WordCloud
from st_aggrid import AgGrid, GridOptionsBuilder
from data.data_fetcher import buscar_artigos, buscar_patentes
from data.database import salvar_interacao

# Configuração inicial
st.set_page_config(page_title="Jornada de Inovação em Saúde", layout="wide")
st.markdown("""
<style>
h1 {color: #1f77b4; font-family: 'Segoe UI', sans-serif;}
h2 {color: #4e89ae;}
.metric-label {font-size: 18px; font-weight: 600;}
.stButton>button {background-color: #1f77b4; color:white; border-radius:5px;}
.expanderHeader {background-color: #e6f0fa; border-radius:5px;}
</style>
""", unsafe_allow_html=True)

st.markdown("<h1 style='text-align: center;'>Jornada Interativa de Inovação</h1>", unsafe_allow_html=True)
st.markdown("<p style='text-align: center; color: #555;'>Explore oportunidades de inovação em saúde de forma rápida e visual.</p>", unsafe_allow_html=True)
st.markdown("---")

# Barra lateral
with st.sidebar:
    st.header("Configurações da Jornada")
    tema = st.text_input("Digite o tema de saúde", placeholder="Ex.: diabetes, vacinas, insulina")
    foco = st.selectbox("Foco", ["Todos", "Tratamento", "Diagnóstico", "Prevenção"])
    st.subheader("Filtro de Datas")
    col1, col2 = st.columns(2)
    data_inicio = col1.number_input("Ano inicial", min_value=1900, max_value=2025, value=2015)
    data_fim = col2.number_input("Ano final", min_value=1900, max_value=2025, value=2025)
    iniciar = st.button("Iniciar Jornada")

# Buscar dados
if iniciar and tema.strip():
    foco_param = None if foco == "Todos" else foco
    df_artigos = buscar_artigos(tema, max_itens=15, foco=foco_param, data_inicio=data_inicio, data_fim=data_fim)
    df_patentes, gaps, recomendacoes = buscar_patentes(tema, df_artigos, max_itens=5)

    st.session_state.df_artigos = df_artigos
    st.session_state.df_patentes = df_patentes
    st.session_state.gaps = gaps
    st.session_state.recomendacoes = recomendacoes
    st.session_state.tema = tema
    st.session_state.foco = foco
    st.session_state.data_inicio = data_inicio
    st.session_state.data_fim = data_fim

    salvar_interacao(f"{tema} - {foco}", len(df_artigos), len(df_patentes))

# Conteúdo principal
if 'df_artigos' in st.session_state:
    df_artigos = st.session_state.df_artigos
    df_patentes = st.session_state.df_patentes
    gaps = st.session_state.gaps
    recomendacoes = st.session_state.recomendacoes
    tema = st.session_state.tema
    foco = st.session_state.foco
    data_inicio = st.session_state.data_inicio
    data_fim = st.session_state.data_fim

    tabs = st.tabs(["Dashboard", "Detalhes", "Relatório"])

    # Dashboard
    with tabs[0]:
        col1, col2, col3 = st.columns(3)
        col1.metric("Artigos encontrados", len(df_artigos))
        col2.metric("Patentes encontradas", len(df_patentes))
        col3.metric("Tema pesquisado", tema)
        st.markdown(f"Período pesquisado: **{data_inicio} - {data_fim}**")
        st.markdown("---")

        grafico_col, nuvem_col = st.columns([1,1])
        with grafico_col:
            st.markdown("Gráfico comparativo")
            tipo_grafico = st.selectbox("Escolha o tipo de gráfico", ["Barras Verticais", "Barras Horizontais"], key="grafico_interativo")
            dados = pd.DataFrame({"Categoria":["Artigos","Patentes"],"Quantidade":[len(df_artigos),len(df_patentes)]})
            if tipo_grafico == "Barras Verticais":
                fig = px.bar(dados, x="Categoria", y="Quantidade", color="Categoria", text="Quantidade",
                             color_discrete_sequence=["#4e89ae","#c56183"])
            else:
                fig = px.bar(dados, y="Categoria", x="Quantidade", orientation="h", color="Categoria", text="Quantidade",
                             color_discrete_sequence=["#4e89ae","#c56183"])
            fig.update_traces(textposition="auto")
            st.plotly_chart(fig, use_container_width=True)

        with nuvem_col:
            st.markdown("Palavras-chave mais frequentes")
            if not df_artigos.empty:
                texto_artigos = " ".join(df_artigos['Título'].tolist())
                wordcloud = WordCloud(width=400, height=300, background_color='#f9f9f9', colormap='Blues').generate(texto_artigos)
                st.image(wordcloud.to_array())
            else:
                st.warning("Nenhum artigo para gerar a nuvem de palavras.")

        st.markdown("---")
        st.subheader("Análise de Oportunidade")
        if len(df_artigos) > 5 and len(df_patentes) < 5:
            st.success("Alta oportunidade de inovação")
        elif len(df_patentes) > 5:
            st.info("Mercado competitivo")
        else:
            st.warning("Tema pouco explorado ou com poucos dados disponíveis")

        # Gaps de mercado
        st.markdown("---")
        st.subheader("Gaps de Mercado")
        if gaps:
            dados_gaps = pd.DataFrame({"Termos": gaps, "Origem": ["Artigos (forte)"] * len(gaps)})
            fig_gaps = px.bar(dados_gaps, x="Termos", y=[1]*len(dados_gaps),
                              color="Origem", text="Origem", title="Principais Gaps Identificados")
            fig_gaps.update_yaxes(visible=False, showticklabels=False)
            st.plotly_chart(fig_gaps, use_container_width=True)

            # Heatmap
            if not df_artigos.empty:
                termo_artigos = pd.Series(" ".join(df_artigos['Título']).lower().split())
                termo_patentes = pd.Series(" ".join(df_patentes['Título']).lower().split()) if not df_patentes.empty else pd.Series([])
                termos_unicos = termo_artigos.unique().tolist()
                termos_comuns = termos_unicos[:15]  # limitar
                dados_heatmap = pd.DataFrame({
                    "Termo": termos_comuns,
                    "Artigos": [sum(termo_artigos == t) for t in termos_comuns],
                    "Patentes": [sum(termo_patentes == t) for t in termos_comuns] if not termo_patentes.empty else [0]*len(termos_comuns)
                })
                fig_heat = px.imshow(dados_heatmap.set_index("Termo").T,
                                     color_continuous_scale="Blues",
                                     title="Heatmap Artigos x Patentes")
                st.plotly_chart(fig_heat, use_container_width=True)

            st.markdown("### Recomendações Estratégicas")
            for rec in recomendacoes:
                st.markdown(f"- {rec}")
        else:
            st.info("Nenhum gap relevante identificado.")

    # Detalhes
    with tabs[1]:
        with st.expander("Detalhes dos Artigos", expanded=True):
            if not df_artigos.empty:
                gb = GridOptionsBuilder.from_dataframe(df_artigos)
                gb.configure_pagination(paginationAutoPageSize=True)
                gb.configure_default_column(sortable=True, filter=True, resizable=True)
                AgGrid(df_artigos, gridOptions=gb.build(), height=250)
            else:
                st.warning("Nenhum artigo encontrado.")

        with st.expander("Detalhes das Patentes", expanded=True):
            if not df_patentes.empty:
                gb = GridOptionsBuilder.from_dataframe(df_patentes)
                gb.configure_pagination(paginationAutoPageSize=True)
                gb.configure_default_column(sortable=True, filter=True, resizable=True)
                AgGrid(df_patentes, gridOptions=gb.build(), height=250)
            else:
                st.warning("Nenhuma patente encontrada.")

    # Relatório 
    with tabs[2]:
        relatorio = f"""
Tema: {tema.title()}
Foco: {foco}
Período: {data_inicio} - {data_fim}
Artigos encontrados: {len(df_artigos)}
Patentes encontradas: {len(df_patentes)}
Análise de Oportunidade: {"Alta oportunidade" if len(df_artigos) > 5 and len(df_patentes) < 5 else "Mercado competitivo" if len(df_patentes) > 5 else "Poucos dados"}

Gaps de Mercado:
{", ".join(gaps) if gaps else "Nenhum gap identificado"}

Recomendações:
{"; ".join(recomendacoes)}

Gerado em: {pd.Timestamp.now().strftime('%d/%m/%Y %H:%M')}
"""
        st.download_button("Baixar relatório (TXT)", data=relatorio, file_name=f"relatorio_{tema}.txt")
        st.markdown("Relatório inclui análise de gaps de mercado e recomendações estratégicas.")

    # Reset 
    if st.button("Nova Jornada"):
        st.session_state.clear()

