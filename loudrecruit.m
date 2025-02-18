
function y = loudrecruit(x)
X = hilbert(x);
env = abs(X);
new_env = env*1.9;
tfs = cos(angle(X));
y = new_env.*tfs;
end