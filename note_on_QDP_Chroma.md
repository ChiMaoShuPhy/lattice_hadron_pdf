
# Link Variables

The relations between link variables stored in Chroma ( $\mathcal{U}$ ) and in text books ( $U$ ) are

 $$ 
 U_\mu(n) = \mathcal{U}_{\mu}(n)
 $$ 



# QDP++ Shift

 The QDP++ shift acts as

 $$ 
 G'=\text{shift}(G,\pm1,\mu) 
 \Rightarrow G'[\mu][n]= G[\mu][n\pm\hat \mu] 
 $$

# Fermion Propagator $S[n]_{\alpha\beta}^{ab}$

The fermion propagator takes the form

$$S[n]_{\alpha_{n}\beta_{src}}^{a_{n}b_{src}}$$ 
where $\alpha,\,\beta$ are Dirac Spinor indices and $a,\,b$ are color indices. $S[n]_{\alpha\beta}^{ab}$ represents the propagator from a lattice site $n$ propagating  ***from*** source site (an input parameter in inverter function)


The gauge transformation of a fermion propagator is

$$S[n]_{\alpha_{n}\beta_{src}}^{a_{n}b_{src}}\rightarrow \sum_{a'=0,b'=0}^{3}g(n)^{a,a'} S[n]_{\alpha_{n}\beta_{src}}^{a'_{n} b'_{src}} g^{\dagger}(src)^{b'b}$$ 

where $g$ is the SU(3) transformation matrix.

S[n] is stored in memory as

$$
S(\{n_x,n_y,n_z,n_t\})\rightarrow\left(
\begin{array}{cc} 
\left(
\begin{array}{cc} 
T_{00}
\end{array}
\right)^{a_{n} b_{src}} & \cdots & \cdots & \left(
\begin{array}{cc} 
T_{03}
\end{array}
\right)^{a_{n} b_{src}} \\
\cdots&\cdots&\cdots&\cdots\\
\cdots&\cdots&\cdots&\cdots\\
\left(
\begin{array}{cc} 
T_{30}
\end{array}
\right)^{a_{n} b_{src}} &\cdots & \cdots &\left(
\begin{array}{cc} 
T_{33}
\end{array}
\right)^{a_{n} b_{src}}
\end{array}
\right)_{\alpha_n\beta_{src}}
$$

where $T_{\alpha_n \beta_{src}}$ is a $3\times 3$ color matrix.

# F-H Source

$$
\left(S(n)_{\alpha'_n \beta'_{src}}^{a'_n b'_{src}}\right)^\dagger \mathcal{L}_{\mu}[n,n+\Delta\hat\mu]_{\alpha_n,\alpha'_{n+\Delta\hat\mu}}^{a_n b_{n+\Delta\hat\mu}}
$$

