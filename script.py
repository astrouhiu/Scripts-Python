# Para importar funções do módulo fibonacci no código atual, escrevemos

import fibonacci

# fib e fib2 são duas funções diferentes que calculam a sequência de Fibonacci

fibonacci.fib(100) 

fibonacci.fib2(1000)

# Outra alternativa à sintaxe de acima é atribuir uma função para uma variável
# com Python

meu_fibonacci = fibonacci.fib

# e a partir daqui usamos simplesmente meu_fibonacci

meu_fibonacci(1000)
