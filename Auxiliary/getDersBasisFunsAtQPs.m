function QData = getDersBasisFunsAtQPs(p, n, KntVect, d, LAB)

uqKntVect = KntVect([true, diff(KntVect) > 0]);
switch LAB
    case 'Gaussian'
        NGPs = (p + 1);
        [Xg, Wg] = GaussRule(NGPs);
        J2 = diff(uqKntVect) / 2;
        QPs = bsxfun(@plus, bsxfun(@times, (Xg + 1), J2'), uqKntVect(1 : end - 1)');
        QPs = reshape(QPs', [], 1);
        W = bsxfun(@times, J2', Wg);
        W = reshape(W', [], 1);
        span = p + 1 : n;
        span = span(diff(KntVect(p + 1 : n + 1)) > 1e-12);
        span = repmat(span, NGPs, 1);
        span = uint32(span(:));
    case 'Optimal'
        uqKntsIdcs = diff(KntVect) > 0;
        KntMult =  diff(find(uqKntsIdcs == 1));
        assert(all(KntMult == KntMult(1)), 'Internal regularity is not consistent, can not use optimal quadrature!')
        r = p - KntMult(1);
        % p = 2, 3, 4, 5
        % even degree: (2p, p ? 2); odd degree: (2p ? 1, p ? 2)
        % Even degree rules provide full integration of isogeometric mass and
        % stiffness matrices, with significantly fewer quadrature points, up to
        % a factor of five in 3D, as compared with Gauss-Legendre quadrature. The
        % odd degree rules can be effectively used in the context of reduced
        % integration, providing another factor of two in computational savings.
        assert( ((p == 2) || (p == 3) || (p == 4) || (p == 5)), 'The degree %d is not supported by optimal integration yet!', p )
        assert( ((p - 2) == (r - 1)), 'The regularity is not supported by optimal integration yet! (Only C^(p - 1) basis functions are currently supported, i.e. maximum regularity)' )
        
        % even degree: full integration of isogeometric mass and stiffness matrices
        [QPs, W] = quadrule(2*p, p - 2, numel(uqKntVect) - 1, KntVect(1), KntVect(end));
        
        %[QPs, W] = quadrule(2*p, r - 1, numel(uqKntVect) - 1, KntVect(1), KntVect(end));
        span = uint32(FindSpan(n, p, QPs, KntVect));
    otherwise
        error('Not implemented yet!');
end

N = mDersBasisFuns(span, QPs, p, d, KntVect); % [p + 1, d + 1, NGPs]

tmp = [true; diff(span) > 0];
spanMult = diff([find(tmp == 1); numel(span) + 1]);
QPIdcs = cumsum([1; spanMult]); % Quadrature point's indices

ElIdcs = zeros(size(QPIdcs));
for i = 1 : numel(spanMult)
    idcs = QPIdcs(i) : QPIdcs(i + 1) - 1;
    ElIdcs(idcs) = i;
end

QData.QPs = QPs;
QData.NumQPs = numel(QPs);
QData.Span = span;
QData.W = W;
QData.N = N;
QData.QPIdcs = uint32(QPIdcs);
QData.ElIdcs = ElIdcs;
end