module DN1

export RazprsenaMatrika
export sor

import Base.lastindex
import Base.firstindex
import Base.getindex
import Base.setindex!
import Base.*

using LinearAlgebra

"""
RazprsenaMatrika(V, I)
    
Podatkovni tip RazprsenaMatrika. V matriki 'V' so nenicelni elementi originalne matrike. V matriki 'I' so 
shranjeni indeksi stolpcev teh nenicelnih elementov. Vsaka vrstica matrike 'V' vsebuje nenicelne elemente iz iste vrstice 
kot originalna matrika.
"""

mutable struct RazprsenaMatrika
    V
    I::Matrix{Int64}
end

"""
RazprsenaMatrika(M)

Funkcija, ki vhodno matriko M pretvori v razprseno matriko. 
Vhod je diagonalno dominantna razprsena matrika, izhod sta matriki 'V' in 'I'.
"""

function RazprsenaMatrika(M)
    max = 1
    for i=1:size(M, 1)
        stevec = 0
        for j=1:length(M[i, :])
            if M[i, j] != 0
                stevec = stevec + 1
            end
        end
        if stevec > max
            max = stevec
        end
    end

    stolpci = max
  
    V = zeros(size(M, 1), stolpci)
    I = zeros(size(M, 1), stolpci)
    vr = 1
    for i=1:size(M, 1)
        count = 1
        for j=1:length(M[i, :])
            if M[i, j] != 0
                V[vr, count] = M[i, j]
                I[vr, count] = j
                count = count + 1
            end
        end
        vr = vr + 1
    end       

    return RazprsenaMatrika(V, I)
end

function getindex(A::RazprsenaMatrika, i::Int64, j::Int64)
    if i > size(A.V, 1) || j > size(A.V, 1) || i < 1 || j < 1
        throw(BoundsError(A, (i, j)))
    end
    for k=1:size(A.I, 2)
        if A.I[i, k] == j
            return A.V[i, k]
        end
    end

    return 0
end

function firstindex(A::RazprsenaMatrika)
    return 1, 1
end


function lastindex(A::RazprsenaMatrika)
    return size(A.V, 1), size(A.V, 1)
end
 
"""
setindex!(A::RazprsenaMatrika, vrednost, i::Int, j::Int)

Funkcija, ki shrani vrednost na doloceni poziciji v matriki. 
Vhod: matrika tipa RazprsenaMatrika, vrednost, indeks vrstice in indeks stolpca. 

"""

function setindex!(A::RazprsenaMatrika, val, i::Int, j::Int)
    if i > size(A.V, 1) || j > size(A.V, 1) || i < 1 || j < 1
        throw(BoundsError(A, (i, j)))
    end

    for k=1:size(A.I, 2)
        if A.I[i, k] == j
            A.V[i, k] = val
            return A
        end
    end

    if j < A.I[i, 1]
        tempV = zeros(size(A.V, 1), size(A.V, 2) + 1)
        tempI = zeros(size(A.V, 1), size(A.V, 2) + 1)

        for k = 1:size(tempV, 1)
            for n = 1: size(tempV, 2)
                if k == i && n == 1
                    tempV[k, n] = val
                    tempI[k, n] = j
                elseif k == i
                    tempV[k, n] = A.V[k, n-1]
                    tempI[k, n] = A.I[k, n-1]
                elseif k != i && n == size(tempV, 2)
                    tempV[k, n] = 0
                    tempI[k, n] = 0
                else
                    tempV[k, n] = A.V[k, n]
                    tempI[k, n] = A.I[k, n]
                end
            end
        end

        A.V = tempV
        A.I = tempI

        return A
    end

    if j > A.I[i, size(A.I, 2)] && A.I[i, size(A.I, 2)] != 0
        tempV = zeros(size(A.V, 1), size(A.V, 2) + 1)
        tempI = zeros(size(A.V, 1), size(A.V, 2) + 1)

        for k = 1:size(tempV, 1)
            for n = 1: size(tempV, 2)
                if k == i && n == size(tempV, 2)
                    tempV[k, n] = val
                    tempI[k, n] = j
                elseif k != i && n == size(tempV, 2)
                    tempV[k, n] = 0
                    tempI[k, n] = 0
                else
                    tempV[k, n] = A.V[k, n]
                    tempI[k, n] = A.I[n, n]
                end
            end
        end
        A.V = tempV
        A.I = tempI
        return A
    end 

    if A.I[i, 1] == 0
        A.I[i, 1] = j
        A.V[i, 1] = val
    else
        for k=1:size(A.I, 2)
            if A.I[i, k] > j
                A.I[i, k+1] = A.I[i, k]
                A.I[i, k] = j
                A.V[i, k+1] = A.V[i, k]
                A.V[i, k] = val
                break
            end
        end
    end
    return A

end


"""
*(A::RazprsenaMatrika, x::Vector)

Funkcija, ki pomnozi matriko tipa RazprsenaMatrika z vektorjem.
"""

function *(A::RazprsenaMatrika, x::Vector)
    if size(A.V, 1) != length(x) 
        throw(error("Dimenzije se ne ujemajo!"))
    end

    b = zeros(length(x))
    print(b)
    for i=1:size(A.V, 1)
        for j=1:size(A.V, 2)
            if A.I[i, j] != 0
                b[i] = b[i] + A.V[i, j]*x[A.I[i, j]]
            end
        end
    end
    return b
end



"""
sor(A::RazprsenaMatrika, b::Vector, x0::Vector, omega, tol=1e-10)

SOR iteracija za resevanje Ax=b

A - matrika tipa RazprsenaMatrika
b - vektor ki je resitev sistema
x0 - zacetni priblizek
omega - relaksacijski parameter
tol - toleranca (pogoj za ustavitev iteracije)
"""

function sor(A::RazprsenaMatrika, b::Vector, x0::Vector, omega, tol=1e-10)
    max = 1000
    for it=1:max
        for i=1:size(A.V, 1)
            sum = 0
            for k=1:size(A.I, 2)
                if i != A.I[i, k] && A.I[i, k] != 0
                    sum = sum + A.V[i, k]*x0[A.I[i, k]]
                end
            end
            x0[i] = (1 - omega) * x0[i] + omega / getindex(A, i, i) * (b[i] - sum)
        end
        
        if norm(A*x0 - b) < tol
            return x0, it
        end
    end

    return x0, max
    
end

end # module DN1