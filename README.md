
# 🧬 Distância Fuzzy entre Sequências de DNA

Este projeto implementa um cálculo de **distância fuzzy** entre sequências de DNA com base em critérios biológicos e métodos de lógica fuzzy. As distâncias são calculadas sem a necessidade de alinhamento entre sequências, utilizando a **Integral de Sugeno**.

---

## 📥 Entrada

- Arquivo FASTA com as sequências de DNA.

---

## ⚙️ Algoritmo

1. **Cálculo de 4 características fuzzy por sequência**:
   - **GC Content**: proporção de bases G e C.
   - **Entropia de Shannon**: medida da diversidade das bases (normalizada).
   - **Padrão CG**: frequência do dígrafo "CG" na sequência.
   - **Similaridade de k-mers via MinHash**: aproxima a similaridade de Jaccard entre pares de sequências com base nos seus k-mers.

2. **Cálculo da função `h`** entre pares de sequências:
   - \[
     h_i = 1 - |\mu_{1i} - \mu_{2i}|
     \]
   - Mede a similaridade entre os valores fuzzy de cada critério.

3. **Cálculo da medida fuzzy λ (lambda)**:
   - λ é a solução da equação da medida fuzzy de Sugeno:
     \[
     \prod_{i}(1 + \lambda \cdot \mu_i) = 1 + \lambda
     \]
   - A solução é obtida numericamente com `fsolve`.

4. **Cálculo da medida fuzzy de subconjuntos (`μ(S)`)**:
   - Usando a fórmula recursiva da medida de Sugeno.

5. **Cálculo da Integral de Sugeno**:
   - Combina os valores de `h_i` com as medidas fuzzy dos subconjuntos para obter a similaridade final.

6. **Cálculo da distância fuzzy**:
   - \[
     \text{distância} = 1 - \text{Integral de Sugeno}
     \]

---

## 📤 Saídas

- `fuzzy_features.csv`: tabela com os valores fuzzy de cada critério para cada sequência.
- `output.csv`: matriz de distâncias fuzzy entre todas as sequências.

---

## 🔧 Parâmetros Importantes

- **Tamanho do k-mer**:
  ```python
  def kmers(seq, k=2)
