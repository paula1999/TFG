# Criptoanálisis del criptosistema de McEliece clásico mediante algoritmos genéticos

> Trabajo de Fin de Grado. Doble Grado en Ingeniería Informática y Matemáticas. Universidad de Granada

Autora: Paula Villanueva Núñez.

Tutor: Gabriel Navarro Garulo.

## :clipboard: Descripción

El objetivo principal de este trabajo se basa en emplear algoritmos genéticos para realizar un análisis del criptosistema de McEliece clásico. En consecuencia, haremos uso de las técnicas y áreas de las matemáticas, tales como los cuerpos finitos, los anillos de polinomios sobre cuerpos finitos y la complejidad de los problemas y de los algoritmos. En cuanto a las herramientas informáticas, hablaremos de las metaheurísticas basadas en poblaciones, en concreto los algoritmos genéticos. A lo largo del desarrollo de este trabajo, ilustraremos en los ejemplos las implementaciones realizadas en el sistema algebraico computacional SageMath.

## :memo: Contenido

El contenido de este trabajo comienza con un capítulo dedicado a describir y desarrollar las herramientas matemáticas e informáticas que necesitaremos para facilitar el seguimiento del posterior desarrollo. Estudiaremos los conceptos básicos relacionados con la teoría de códigos lineales, que nos permitirá introducir las bases para que en el siguiente capítulo podamos desarrollar los códigos de Goppa. A partir de estos códigos, podremos construir el criptosistema de McEliece clásico y emplear los algoritmos genéticos para concluir con su respectivo análisis.

## :heavy_check_mark: Objetivos

1. Estudiar la teoría básica de códigos lineales.
2. Estudiar los códigos de Goppa y su decodificación.
3. Estudiar la criptografía basada en códigos como modelo de criptografía post-cuántica.
4. Estudiar e implementar el criptosistema de McEliece.
5. Estudiar e implementar algoritmos evolutivos para el cálculo de la distancia de un código lineal.
6. Estudiar el criptoanálisis del criptosistema de McEliece mediante algoritmos evolutivos.

## :book: Memoria

Para generar la memoria se debe descargar el repositorio y usar la herramienta `make` dentro de la carpeta [`documentation`](./documentation).

## :computer: Código

El código implementado se encuentra en la carpeta [`src`](./src). Para usar el código desarrollado primero se debe cargar cada fichero `.sage` con la función `load()` de SageMath.
