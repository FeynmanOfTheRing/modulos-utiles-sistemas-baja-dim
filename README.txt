Aquí hay 3 módulos principales que me parece nos pueden servir


Uno es una clase tipo caja, donde se puede crear un objeto tipo caja cuántica
que tenga como atributos su longitud y su altura de la barrera, teniendo así 
sus autofunciones y sus autovalores de energía.


Otro es una clase tipo oscilador armónico, donde se crea un objeto tipo oscilador
armónico mecánico cuántico que tiene como atributo el factor de parabolicidad hbW, 
de donde también obtenemos los valores de la energía y las autofunciones.


Y el que me parece más útil es el módulo donde están programadas las fórmulas de 
absorción pero con la facilidad de que se pueden ver de manera simbólica como están 
escritas (función: mostrar_formulas()), y así asegurarnos que están bien programadas.
Trabajar con esto nos permite aplicar esta fórmula a cualquier otro problema 
solamente definiendo nuestros parámetros y calculando la diferencia de 
energía entre los niveles 1 y 2 (E21) y las integrales del momento dipolar
(M22,M11,M12), haciendo que el código sea mucho más corto y fácil de usar.