function randnNEW(a::Int64, b::Int64)
    global c
    global mydraws
    d = deepcopy(c + a * b)
    bla = mydraws[c + 1:d]
    c = deepcopy(d)
    return reshape(bla, (a, b))
end