function P = probAtleastOne(p);

P = 1 - cumprod(1-p);