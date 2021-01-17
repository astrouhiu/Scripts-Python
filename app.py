'''
Data App - Previsão do Valor do Imóvel
Para executar digite: streamlit run app.py
'''

import streamlit as st 
import pandas as pd
import plotly.express as px
from sklearn.ensemble import RandomForestRegressor

# Função para carregar o dataset
@st.cache
def get_data():
    return pd.read_csv("data.csv")

# Função para treinar o modelo
def train_model():
    data = get_data()
    # Separa os dados
    x = data.drop("MEDV", axis = 1)
    y = data["MEDV"]
    # Instancia o RandomForestRegressor passando alguns parâmetros
    rf_regressor = RandomForestRegressor(n_estimators = 200, max_depth = 7, max_features = 3)
    # Treina o modelo
    rf_regressor.fit(x, y)
    # Retorna o classificador
    return rf_regressor

# Cria um DataFrame
data = get_data()

# Treina o modelo
model = train_model()

# Título
st.title("Data App - Previsão do Valor do Imóvel")

# Subtítulo
st.markdown("Este Data App exibe a solução de Machine Learning na predição de valores de imóveis da cidade de Boston")

# Verifica o dataset
st.subheader("Mostrando apenas um pequeno conjunto dos atributos")

# Atributos para serem exibidos por padrão
defaultcols = ["RM", "PTRATIO", "LSTAT", "MEDV"]

# Define atributos a partir do multiselect e exibe na tela aqueles em 'default'
cols = st.multiselect("Atributos", data.columns.tolist(), default = defaultcols)

# Exibe os top 10 registros do DataFrame
st.dataframe(data[cols].head(10))

st.subheader("Distribuição de imóveis por preço na escala de US$ 1000")

# Define a faixa de valores
faixa_valores = st.slider("Faixa de preço", float(data.MEDV.min()), 50., (10., 40.))

# Filtra os dados
dados = data[data['MEDV'].between(left = faixa_valores[0], right = faixa_valores[1])]

# Plota a distribuição dos dados com filtragem dinâmica
f = px.histogram(dados, x = "MEDV", nbins = 50, title = "Distribuição de Preços")
f.update_xaxes(title = "MEDV (US$ 1000)")
f.update_yaxes(title = "Total Imóveis")
st.plotly_chart(f)

st.sidebar.subheader("Defina os atributos do imóvel para predição")

# Mapea dados do usuário para cada atributo
crim = st.sidebar.number_input("Taxa de criminalidade", value = data.CRIM.mean())
indus = st.sidebar.number_input("Proporção de hectares de negócio", value = data.CRIM.mean())
chas = st.sidebar.selectbox("Faz limite com o rio?", ("Sim", "Não"))
# Transforma o dado de entrada em valor binário
chas = 1 if chas == "Sim" else 0
nox = st.sidebar.number_input("Concentração de óxido nítrico", value = data.NOX.mean())
rm = st.sidebar.number_input("Número de quartos", value = 1)
ptratio = st.sidebar.number_input("Índice de alunos para professores", value = data.PTRATIO.mean())
b = st.sidebar.number_input("Proporção de pessoas com descendencia afro-americana", value = data.B.mean())
lstat = st.sidebar.number_input("Porcentagem de status baixo", value = data.LSTAT.mean())

# Insere um botão na tela
btn_predict = st.sidebar.button("Realizar predição")

# Verifica se o botão for acionado
if btn_predict:
    result = model.predict([[crim, indus, chas, nox, rm, ptratio, b, lstat]])
    st.subheader("O valor previsto para o imóvel é:")
    result = "US $ " + str(round(result[0]*1000, 2))
    st.write(result)
