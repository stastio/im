function f(x, dyn::Model)
[x[2], -9.8sin(x[1])]
end

function ∇f(x, dyn::Model)
[0 1; -9.8cos(x[1]) 0]
end

function g(x, dyn::Model)
[0, 1]
end

function ∇g(x, dyn::Model)
[0 0; 0 0]
end

function h(x, dyn::Model)
[0, 1]
end

function ∇h(x, dyn::Model)
[0 0; 0 0]
end
