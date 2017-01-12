### Es ist ein MPP "C++"-Projekt...

Mit dem Aufruf:

        mpirun -np 4 project_c -f data8 -g 8

wird eine leere Matrix in die Datei _**data8**_ mit 8 * 8 elementen generiert.

Der Aufruf:

        mpirun -np 4 project_c -f data8 -g 8_2_3_4_5

generiert eine Matrix mit 8 * 8 Elementen mit einem Rechteck in der
Position x=2, y=3, der Länge w=4 und der Höhe h=5.

Die Matrix kann mit dem Flag __-p__ ausgegeben werden.

        mpirun -np 4 project_c -f data8 -p

Ausgabe:

        ........
        ........
        ........
        ..XXXX..
        ..XXXX..
        ..XXXX..
        ..XXXX..
        ..XXXX..

#### Indizierung fängt bei 0 an.

----

n der Größe 10000 bedeutet 10000² + 10000 (_'\\n'_) = **100.010.000 Zeichen**.
Da ein Zeichen 1 Byte beträgt, resultiert sich eine Größe von **100,01 MB**.
100,01 MB entspricht: 100.010.000 / 1024 / 1024 = **95,376968384 MiB**.

Die generierte Dateigröße wächst exponentiell (auf _'\\n'_ am Ende der Zeile
verzichtet).

n           |Größe in Bytes     |Größe
------------|-------------------|-----
1           |1                  |1 B
10          |100                |100 B
100         |10.000             |10 kB
1.000       |1.000.000          |1 MB
10.000      |100.000.000        |**100 MB**
20.000      |400.000.000        |**400 MB**
30.000      |900.000.000        |**900 MB**
40.000      |1.600.000.000      |**1,6 GB**
50.000      |2.500.000.000      |**2,5 GB**
60.000      |3.600.000.000      |**3,6 GB**
70.000      |4.900.000.000      |**4,9 GB**
80.000      |6.400.000.000      |**6,4 GB**
90.000      |8.100.000.000      |**8,1 GB**
100.000     |10.000.000.000     |**10 GB**
1.000.000   |1.000.000.000.000  |**1 TB**

max n = 46340, da n² = 2147395600 (int_max = 2147483647)