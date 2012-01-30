function inv_hb = inv_get_chromophore(h1)
switch h1
    case 'HbO'
        inv_hb = 1;
    case 'HbR'
        inv_hb = 2;
    case 'HbT'
        inv_hb = 3;
end