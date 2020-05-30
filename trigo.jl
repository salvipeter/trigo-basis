module Trigo

"""
    generate_odd(n; pretty)

Generates the coefficients of the odd triangle until the `n`th row.
Print human-readable output when `pretty` is `true`.
"""
function generate_odd(n; pretty=false)
    coeffs = Dict((0,0,0)=>1)
    multiplier = Dict((1,0,0)=>1, (0,1,0)=>1, (0,0,1)=>1)
    result = []
    for m in 1:n
        coeffs = sort(poly_multiply(coeffs, multiplier), rev=true)
        push!(result, map(c -> Dict(c[1]=>c[2]), collect(coeffs)))
        if pretty
            println("\nm = $m:")
            for poly in result[end]
                print("\t")
                poly_print(poly)
            end
        end
    end
    result
end

"""
    generate_even(n; pretty)

Generates the coefficients of the even triangle until the `n`th row.
Print human-readable output when `pretty` is `true`.
"""
function generate_even(n; pretty=false)
    odd_coeffs = generate_odd(n)
    a_hat = Dict((0,0,0)=>1, (0,0,1)=>-2, (0,0,2)=>1)
    c_hat = Dict((0,0,0)=>1, (1,0,0)=>-2, (2,0,0)=>1)
    result = []
    for m in 1:n
        row = odd_coeffs[m]
        coeffs = [row[1]]
        for i in 1:m
            push!(coeffs, poly_add(row[2i], row[2i+1]))
        end
        for i in 1:m
            push!(coeffs, poly_add(row[2i-1], row[2i]))
        end
        push!(coeffs, row[end])
        for i in 1:m+1
            coeffs[i] = poly_multiply(coeffs[i], a_hat)
        end
        for i in m+2:2m+2
            coeffs[i] = poly_multiply(coeffs[i], c_hat)
        end
        push!(result, coeffs)
        if pretty
            println("\nm = $m:")
            for poly in result[end]
                print("\t")
                poly_print(poly)
            end
        end
    end
    result
end

"""
    write_table(n, filename)

Write a computer-readable representation of the triangle into a file.
The odd and even lines are interleaved, and each polynomial is followed
by its first and second derivatives.
"""
function write_table(n, filename)
    odd_rows = generate_odd(n)
    even_rows = generate_even(n)
    open(filename, "w") do f
        println(f, 2 * n)         # number of rows we are going to write
        for m in 1:n
            println(f, 2 * m + 1) # number of terms in the current row
            for term in odd_rows[m]
                poly_write(term, f)
                d = poly_derivative(term)
                poly_write(d, f)
                poly_write(poly_derivative(d), f)
            end
            println(f, 2 * m + 2) # number of terms in the current row
            for term in even_rows[m]
                poly_write(term, f)
                d = poly_derivative(term)
                poly_write(d, f)
                poly_write(poly_derivative(d), f)
            end
        end
    end
end

"""
    poly_add(x, y)

Add two trigonometric polynomials
"""
function poly_add(x, y)
    result = copy(x)
    for (abc, k) in y
        result[abc] = get(result, abc, 0) + k
    end
    result
end

"""
    poly_multiply(x, y)

Multiply two trigonometric polynomials.
"""
function poly_multiply(x, y)
    result = Dict()
    for (abc1, k1) in x, (abc2, k2) in y
        a = abc1[1] + abc2[1]
        b = abc1[2] + abc2[2]
        c = abc1[3] + abc2[3]
        k = k1 * k2
        if b == 2
            a += 1
            b = 0
            c += 1
            k *= 2
        end
        result[(a,b,c)] = get(result, (a,b,c), 0) + k
    end
    result
end

"""
    poly_scale(p, x)

Scales the trigonometric polynomial `p` by the real number `x`.
"""
function poly_scale(p, x)
    result = Dict()
    for (abc, k) in p
        result[abc] = k * x
    end
    result
end

"""
    poly_print(x)

Prints a trigonometric polynomial in human_readable form.
"""
function poly_print(x)
    first = true
    for (abc, k) in x
        !first && k >= 0 && print("+")
        k != 1 && print(k)
        for i in 1:3
            if abc[i] > 0
                print("abc"[i])
                abc[i] > 1 && print("^$(abc[i])")
            end
        end
        first = false
    end
end

"""
    poly_write(p, io::IO)

Writes the trigonometric polynomial `p` onto the I/O stream `io`
in a computer-readable format.
"""
function poly_write(p, io::IO = stdout)
    print(io, length(p))
    for (abc, k) in p
        print(io, " $(abc[1]) $(abc[2]) $(abc[3]) $k")
    end
    println(io)
end

"""
    poly_derivative(x)

Computes the derivative of a trigonometric polynomial.
"""
function poly_derivative(x)
    result = Dict()
    for (abc, k) in x
        result = poly_add(result, derivative_term(abc, k))
    end
    result
end

"""
    derivative_term(abc, k)

Computes the derivative of a single trigonometric term.
The result is a trigonometric polynomial.
"""
function derivative_term(abc, k)
    poly_scale(poly_add(derivative_term(abc, 1, 1),
                        poly_add(derivative_term(abc, 1, 2),
                                 derivative_term(abc, 1, 3))),
               k)
end

"""
    derivative_term(abc, k, i)

Computes the derivative based on just the `i`th subterm.
"""
function derivative_term(abc, k, i)
    abc[i] == 0 && return Dict()
    base = derivative_base(i)
    index = collect(abc)
    index[i] -= 1
    index = (index[1], index[2], index[3])
    poly_multiply(Dict(index=>k*abc[i]), base)
end

"""
    derivative_base(i)

Returns the derivative of the `i`th base.
"""
function derivative_base(i)
    if i == 1
        Dict((0,0,0)=>-1, (0,0,1)=>1) # a' = c - 1
    elseif i == 2
        Dict((1,0,0)=>1, (0,0,1)=>-1) # b' = a - c
    else # i == 3
        Dict((0,0,0)=>1, (1,0,0)=>-1) # c' = 1 - a
    end
end

end # module
