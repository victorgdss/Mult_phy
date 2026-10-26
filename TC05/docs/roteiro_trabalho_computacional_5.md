# MODELAGEM E SIMULAÇÃO MULTIFÍSICA
## TRABALHO COMPUTACIONAL 5
### O levitador magnético — Modelagem usando PDEs

**Objetivo:** Nesta aula, retomaremos o estudo do problema do levitador magnético, agora analisando-o sob a perspectiva das equações de Maxwell. Utilizaremos o método de elementos finitos (FEM) para modelar as grandezas magnéticas, obtendo os valores concentrados necessários para a simulação do modelo dinâmico.

**Exemplos de problemas magnetostáticos:** Para uma melhor compreensão da utilização do FEniCS em problemas magnéticos axissimétricos, siga os tutoriais apresentados nas aulas 13 e 18. O tutorial da aula 13 apresenta uma bobina no espaço livre. Nesse tutorial é apresentada a forma fraca do problema magnetostático axissimétrico e como o resultado pode ser usado para determinar a impedância da bobina. Na aula 18, o problema apresentado na aula 13 é modificado para receber o núcleo ferromagnético e uma peça cilíndrica em levitação. No caso da permeabilidade dos materiais, tanto a permeabilidade magnética do cobre (bobina) quanto a do ar são aproximadas pela permeabilidade magnética do vácuo, $\mu_{0} = 4\pi \cdot 10^{-7} \, (\text{H/m})$. Para as partes ferromagnéticas, são considerados materiais lineares e isotrópicos com $\mu_{r} = 1200$. Reparem que ambos os exemplos são regidos pelas mesmas equações, o que muda é apenas o domínio de interesse.

Ainda na aula 18, é apresentado um levitador de geometria mais simples, cuja solução analítica torna-se viável ao considerar infinita a permeabilidade das partes ferromagnéticas. A partir da solução desse problema, podem ser determinadas as grandezas $L$, $Bl$ e $f_m$, utilizadas na simulação do modelo dinâmico. Como esperado, esses parâmetros dependem da posição da peça em levitação, da corrente e de uma constante $k$, que é obtida a partir dos parâmetros geométricos do problema.

---

### Geometria do problema

Considere a geometria apresentada na Fig. 1, onde os parâmetros construtivos do levitador estão descritos na Tabela 1.

| Descrição | Parâmetro | Valor |
| :--- | :---: | :--- |
| Raio do núcleo | $r_c$ | $0.75 \, \text{cm}$ |
| Altura do núcleo | $z_c$ | $6.0 \, \text{cm}$ |
| Distância entre o núcleo e a esfera | $T$ | $0.5 \, \text{cm}$ |
| Raio da esfera | $r_s$ | $1.0 \, \text{cm}$ |
| Massa da esfera | $m$ | $22.0 \, \text{g}$ |
| Altura da bobina | $z_h$ | $4.0 \, \text{cm}$ |
| Bitola do fio | AWG | $18$ |
| Número de espiras | $N$ | $200$ |

*TABLE 1. Parâmetros do levitador*

#### Descrição da figura 1:

Com base no documento fornecido, a geometria do levitador magnético (ilustrada na Figura 1 do roteiro) consiste num núcleo cilíndrico fixo, uma bobina condutora enrolada à sua volta, e uma esfera ferromagnética situada abaixo do núcleo que se encontra em levitação.   

O problema será resolvido a partir de um modelo modelo 2D (Axissimétrico), que simplifica o domínio do problema para um plano $r$-$z$ explorando a simetria cilíndrica, o que otimiza o cálculo computacional através do método dos elementos finitos (FEM). 

O modelo 2D axissimétrico é constituído por um domínio delimitado por um retângulo maior que representa o espaço livre (ar). No canto superior esquerdo, encostado ao eixo de rotação ($r = 0$), encontra-se o retângulo correspondente ao núcleo cilíndrico de ferro. A bobina é representada por um retângulo adjacente à lateral externa do núcleo. Abaixo do conjunto núcleo-bobina, respeitando a distância de entreferro especificada, posiciona-se a semiesfera ferromagnética, com seu lado plano apoiado sobre o eixo de simetria.

---

### Aproximação Linear

Inicialmente, assumiremos que os materiais ferromagnéticos apresentam um comportamento linear. Nesse contexto, consideraremos que tanto o núcleo quanto a esfera possuem $\mu_r = 1200$.

#### Parte prática

Simule o problema considerando os parâmetros descritos na Tab. 1. Inicialmente, utilize uma corrente de $1 \, \text{A}$. A seguir, descreva os passos utilizados na simulação: definição da forma fraca, tipo de condições de contorno, cálculo da densidade de corrente, entre outros.

Apresente os resultados do Campo Magnético ($H$) e da Indução Magnética ($B$) dentro do domínio. Além disso, mostre as curvas das componentes $H_r$ e $H_z$ na superfície externa da esfera (campos avaliados fora do material ferromagnético). A partir dessas componentes e do vetor normal à superfície da esfera, calcule a componente da força magnética por meio do tensor de Maxwell.

Por fim, compare esse resultado com o obtido utilizando as funções `FacetNormal()` e `assemble()` do FEniCS. Ajuste o valor da corrente de forma que a força magnética compense exatamente a força peso. Esse valor é compatível com o tipo de fio utilizado? Qual a potência dissipada no fio?

---

### Levantamento dos Parâmetros de Circuito

Com a corrente corrigida, determine os parâmetros concentrados utilizados no modelo circuital, considerando a posição de equilíbrio e outros pontos linearmente espaçados próximos a essa condição. Utilize como referência os dados apresentados na Tabela 2. Para esses pontos, anote também os valores da Indução Magnética máxima e média dentro do núcleo e dentro da esfera. Organize seus resultados em uma ou mais tabelas.

Levante as curvas de indutância e força magnética em função da posição. Compare os resultados obtidos com as expressões analíticas derivadas para o modelo do levitador planar apresentado na Aula 18.

| Parâmetro | Valor |
| :--- | :--- |
| Ponto inicial | $0.1$ |
| Ponto final | $3.0$ |
| Número de pontos | $10$ |

*TABLE 2. Pontos para a análise dos parâmetros concentrados.*

---

### Efeitos da Não Linearidade

A inclusão dos efeitos da não linearidade no comportamento dos materiais ferromagnéticos altera significativamente os resultados da análise. Em particular, a permeabilidade relativa ($\mu_r$) deixa de ser constante, tornando-se uma função da densidade de fluxo magnético ($B$). Esse comportamento resulta em uma relação não linear entre o campo magnético ($H$) e a indução magnética ($B$), conforme exemplificado na figura 2.

Os dados apresentados na figura 2 foram obtidos da biblioteca de materiais do software FEMM e são apresentados na tabela 3.

No contexto da simulação, a não linearidade implica a necessidade de utilizar métodos numéricos iterativos para a solução do problema, uma vez que a forma fraca resultante deixa de ser linear. Além disso, os valores de indutância e força magnética passam a depender não apenas da posição da peça, mas também da intensidade da corrente aplicada, devido à saturação magnética nas regiões de maior densidade de fluxo.

Dependendo do ponto de equilíbrio e do deslocamento da esfera, esses efeitos devem ser considerados na análise e no projeto do sistema, a fim de garantir que os modelos empregados representem adequadamente o comportamento real do dispositivo.

Com base nos valores de $B_{\text{max}}$ e $B_{\text{med}}$ avaliados anteriormente, determine os pontos de operação correspondentes na curva $B \times H$. Analise a diferença entre o valor esperado da permeabilidade no problema não linear em comparação com a aproximação linear.

---

| Índice | $B$ ($\text{Tesla}$) | $H$ ($\text{A/m}$) |
| :---: | :--- | :--- |
| 1 | $0.000000$ | $0.0000$ |
| 2 | $0.227065$ | $13.8984$ |
| 3 | $0.454130$ | $27.7967$ |
| 4 | $0.681195$ | $42.3974$ |
| 5 | $0.908260$ | $61.4157$ |
| 6 | $1.135330$ | $82.3824$ |
| 7 | $1.362390$ | $144.6690$ |
| 8 | $1.589350$ | $897.7600$ |
| 9 | $1.812360$ | $4581.7400$ |
| 10 | $2.010040$ | $17736.0000$ |
| 11 | $2.133160$ | $41339.3000$ |
| 12 | $2.199990$ | $68321.8000$ |
| 13 | $2.254790$ | $95685.5000$ |
| 14 | $2.299930$ | $123355.0000$ |
| 15 | $2.342510$ | $151083.0000$ |
| 16 | $2.378760$ | $178954.0000$ |
| 17 | $2.415010$ | $206825.0000$ |
| 18 | $2.451260$ | $234696.0000$ |
| 19 | $2.487500$ | $262568.0000$ |
| 20 | $2.523750$ | $290439.0000$ |
| 21 | $2.560000$ | $318310.0000$ |

*TABLE 3. Valores da curva de magnetização do ferro puro.*

---

**Observação de Entrega:** Escreva um relatório contendo as atividades desenvolvidas. Converta o arquivo para o formato PDF e o envie pelo Moodle. Anexe os arquivos de simulação em um único arquivo zip e envie juntamente com o PDF. Cuidado! Relatórios que não estiverem em formato PDF não serão avaliados.