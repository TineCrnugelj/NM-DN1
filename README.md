# SOR iteracija za razpršene matrike

Osnovni tip je RazprsenaMatrika, ki shranjuje neničelne vrednosti matrike v matriki `V` ter indekse teh vrednosti v matriki `I`.
Funkcija `getindex` omogoča dostop do elementa na določenem indeksu, medtem ko `setindex` dodaja nove vrednosti v matriko. Ključna funkcija `sor` (Successive Over-Relaxation) izvaja iteracijo za reševanje sistema linearnih enačb Ax = b, kjer je A razpršena matrika. SOR je iterativna metoda, ki uporablja relaksacijski parameter omega za hitrejšo konvergenco k rešitvi. Iteracija se ustavi, ko je dosežena določena toleranca napake.
