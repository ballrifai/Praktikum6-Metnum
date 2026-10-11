clc; clear; close all;
format long g

A = [ 5 -1 -1  0  0;
     -1  5 -1 -1  0;
     -1 -1  4 -1 -1;
      0  0  1  4 -2;
      0  1 -1  1  4];
b  = [-1; 2; 6; 2; -1];
n  = length(b);
x0 = zeros(n,1);

nRep    = 2000;
nIter   = 2;
tol     = 1e-10;
maxIter = 1000;
nama    = {'p','q','r','s','t'};

%% ====================== 1. ELIMINASI GAUSS ===========================
tic
for rep = 1:nRep
    Ab = [A b];
    for k = 1:n-1
        for i = k+1:n
            m = Ab(i,k) / Ab(k,k);
            Ab(i,k:n+1) = Ab(i,k:n+1) - m*Ab(k,k:n+1);
        end
    end
    xG = zeros(n,1);
    xG(n) = Ab(n,n+1) / Ab(n,n);
    for i = n-1:-1:1
        xG(i) = (Ab(i,n+1) - Ab(i,i+1:n)*xG(i+1:n)) / Ab(i,i);
    end
end
tGauss = toc / nRep;

%% ==================== 2. ELIMINASI GAUSS-JORDAN ======================
tic
for rep = 1:nRep
    Ab = [A b];

    for k = 1:n
        Ab(k,:) = Ab(k,:) / Ab(k,k);

        for i = 1:n
            if i ~= k
                m = Ab(i,k);
                Ab(i,:) = Ab(i,:) - m * Ab(k,:);
            end
        end
    end

    xGJ = Ab(:,n+1);
end
tGJ = toc / nRep;

%% ================== 3. DEKOMPOSISI LU (DOOLITTLE) ====================
tic
for rep = 1:nRep
    L = eye(n);  U = zeros(n);
    for i = 1:n
        for j = i:n
            U(i,j) = A(i,j) - L(i,1:i-1)*U(1:i-1,j);
        end
        for j = i+1:n
            L(j,i) = (A(j,i) - L(j,1:i-1)*U(1:i-1,i)) / U(i,i);
        end
    end
    y = zeros(n,1);
    for i = 1:n
        y(i) = b(i) - L(i,1:i-1)*y(1:i-1);
    end
    xLU = zeros(n,1);
    for i = n:-1:1
        xLU(i) = (y(i) - U(i,i+1:n)*xLU(i+1:n)) / U(i,i);
    end
end
tLU = toc / nRep;

%% ============ 4a. JACOBI - TEPAT 2 ITERASI (sesuai manual) ===========
tic
for rep = 1:nRep
    x = x0;
    for k = 1:nIter
        xn = x;
        for i = 1:n
            idx = [1:i-1, i+1:n];
            xn(i) = (b(i) - A(i,idx)*x(idx)) / A(i,i);
        end
        x = xn;
    end
end
tJ2 = toc / nRep;
xJ2 = x;

%% ========== 4b. JACOBI - SAMPAI KONVERGEN (|dx|_inf < tol) ===========
tic
for rep = 1:nRep
    x = x0;
    for k = 1:maxIter
        xn = x;
        for i = 1:n
            idx = [1:i-1, i+1:n];
            xn(i) = (b(i) - A(i,idx)*x(idx)) / A(i,i);
        end
        selisih = norm(xn - x, inf);
        x = xn;
        if selisih < tol, break; end
    end
end
tJc   = toc / nRep;
xJc   = x;
iterJ = k;

%% ========= 5a. GAUSS-SEIDEL - TEPAT 2 ITERASI (sesuai manual) ========
tic
for rep = 1:nRep
    x = x0;
    for k = 1:nIter
        for i = 1:n
            idx = [1:i-1, i+1:n];
            x(i) = (b(i) - A(i,idx)*x(idx)) / A(i,i);
        end
    end
end
tGS2 = toc / nRep;
xGS2 = x;

%% ======= 5b. GAUSS-SEIDEL - SAMPAI KONVERGEN (|dx|_inf < tol) ========
tic
for rep = 1:nRep
    x = x0;
    for k = 1:maxIter
        xlama = x;
        for i = 1:n
            idx = [1:i-1, i+1:n];
            x(i) = (b(i) - A(i,idx)*x(idx)) / A(i,i);
        end
        if norm(x - xlama, inf) < tol, break; end
    end
end
tGSc   = toc / nRep;
xGSc   = x;
iterGS = k;

xEks = xG;     %

fprintf('\n================ HASIL METODE LANGSUNG ================\n');
fprintf('%-4s %16s %16s %16s\n','Var','Gauss','Gauss-Jordan','LU');
for i = 1:n
    fprintf('%-4s %16.10f %16.10f %16.10f\n', nama{i}, xG(i), xGJ(i), xLU(i));
end

% --- Riwayat iterasi (untuk mencocokkan perhitungan manual) ---
xj = x0;  xg = x0;
histJ = zeros(n,nIter);  histGS = zeros(n,nIter);
for k = 1:nIter
    xn = xj;
    for i = 1:n
        idx = [1:i-1, i+1:n];
        xn(i) = (b(i) - A(i,idx)*xj(idx)) / A(i,i);
    end
    xj = xn;  histJ(:,k) = xj;
    for i = 1:n
        idx = [1:i-1, i+1:n];
        xg(i) = (b(i) - A(i,idx)*xg(idx)) / A(i,i);
    end
    histGS(:,k) = xg;
end
fprintf('\n================ ITERASI MANUAL ================\n');
fprintf('%-4s %12s %12s | %12s %12s\n','Var','Jacobi-1','Jacobi-2','GS-1','GS-2');
for i = 1:n
    fprintf('%-4s %12.8f %12.8f | %12.8f %12.8f\n', nama{i}, ...
        histJ(i,1), histJ(i,2), histGS(i,1), histGS(i,2));
end

% --- Tabel galat mutlak iterasi ke-2 ---
errJ  = abs(xEks - xJ2);
errGS = abs(xEks - xGS2);
fprintf('\n================ GALAT MUTLAK ITERASI KE-2 ================\n');
fprintf('%-4s %13s | %13s %13s | %13s %13s\n', 'Var','Eksak', ...
        'Jacobi-2','Galat','GS-2','Galat');
for i = 1:n
    fprintf('%-4s %13.8f | %13.8f %13.8f | %13.8f %13.8f\n', nama{i}, ...
        xEks(i), xJ2(i), errJ(i), xGS2(i), errGS(i));
end
fprintf('Galat maks (norma tak hingga): Jacobi = %.8f | Gauss-Seidel = %.8f\n', ...
        max(errJ), max(errGS));

% --- Konvergensi sampai toleransi ---
fprintf('\n================ ITERASI SAMPAI KONVERGEN (tol = %g) ================\n', tol);
fprintf('Jacobi       : %d iterasi, galat maks = %.3e\n', iterJ,  max(abs(xEks - xJc)));
fprintf('Gauss-Seidel : %d iterasi, galat maks = %.3e\n', iterGS, max(abs(xEks - xGSc)));

% --- Radius spektral (indikator laju konvergensi) ---
D  = diag(diag(A));  Lo = tril(A,-1);  Up = triu(A,1);
rhoJ  = max(abs(eig(-D \ (Lo + Up))));
rhoGS = max(abs(eig(-(D + Lo) \ Up)));
fprintf('\nRadius spektral: Jacobi = %.6f | Gauss-Seidel = %.6f (harus < 1)\n', rhoJ, rhoGS);

% --- Running time ---
fprintf('\n================ RUNNING TIME (rata-rata dari %d pengulangan) ================\n', nRep);
fprintf('%-34s %14s\n','Metode','Waktu (mikrodetik)');
fprintf('%-34s %14.3f\n','Eliminasi Gauss',                tGauss*1e6);
fprintf('%-34s %14.3f\n','Eliminasi Gauss-Jordan',         tGJ*1e6);
fprintf('%-34s %14.3f\n','Dekomposisi LU',                 tLU*1e6);
fprintf('%-34s %14.3f\n','Jacobi (2 iterasi)',             tJ2*1e6);
fprintf('%-34s %14.3f\n',sprintf('Jacobi (konvergen, %d it.)', iterJ),   tJc*1e6);
fprintf('%-34s %14.3f\n','Gauss-Seidel (2 iterasi)',       tGS2*1e6);
fprintf('%-34s %14.3f\n',sprintf('Gauss-Seidel (konvergen, %d it.)', iterGS), tGSc*1e6);
