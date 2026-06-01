function q = q_product(a, b)
% Prodotto di due quaternioni 4x1, formato scalar-last [qv; qs]
% a, b: [4x1] column vectors
% q:  [4x1] column vector = a * b
av = a(1:3);
as = a(4);
bv = b(1:3);
bs = b(4);

qv = as * bv + bs * av + cross(av, bv);
qs = as * bs - dot(av, bv);

q = [qv; qs];
end