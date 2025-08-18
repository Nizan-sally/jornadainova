# teste_salvar.py
from data.database import salvar_interacao

# Força uma interação de teste
salvar_interacao("teste manual", total_artigos=2, total_patentes=1)