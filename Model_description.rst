===============
Microstructure
===============

Key features of the 316H stainless steel microstructure are noted below:
1. The creep behavior of 316H shows multiple secondary stages, which also appears to be dependent on the microstructural state.
2. Average grain size was observed to be around 50 μm.
3. The primary precipitates/particles found in 316H are Cr23C6 (carbide) and Fe2Mo (Laves).
4. The particles typically grow with deformation.
5. The primary solute particles are C, Cr and Mo, the concentrations of which typically reduce with deformation


Flow Rule
-------------
The deformation kinetics on a given slip system is given by the following expression
:math:`\dot{\gamma} = \dot{\gamma_0}\left|\frac{\tau+\tau_{in}}{\tau_{cr}}\right|^p sgn(\tau+\tau_{in})`

..
:math:`\dot{\gamma} = \dot{\gamma_0} \exp \left( -\frac{\Delta F_0}{kT} \left(1-\left|\frac{\tau + \tau_{in}}{(1+F_{sol})\tau_{cr}}\right|^{3/4}\right)^{4/3} \right) sgn(\tau)`


The most basic ODEs are in given in the forms of :math:`\dot{L_d}`, :math:`\dot{c}`, :math:`\dot{r_p}` and :math:`\dot{f_v}`.

:math:`\dot{L_d} = -L_{d}^3 \left( J_1\dot{\gamma_i} + J_2 \sum_{j \neq i} \dot{\gamma_j}  \right) + K\frac{1}{L_{d}^3}`

:math:`\dot{c} = N_0Z\beta \exp\left(-\frac{G^*}{kT}\right)\frac{4\pi}{3} \left(\frac{D^r}{r_p} \frac{c^r-c^r_{eq}}{c_p^r-c^r_{eq}} + \frac{N_0 Z \beta \exp {\left( - \frac{G^*}{kT} \right)}}{N_v} (r_c - r_p)\Delta t \right)^3 \frac{c_0-c_p}{(1-f_v)}`

:math:`\dot{f_v} = (1-f_v) N_0Z\beta \exp\left(-\frac{G^*}{kT}\right)\frac{4\pi}{3} \left(\frac{D^r}{r_p} \frac{c^r-c^r_{eq}}{c_p^r-c^r_{eq}} + \frac{N_0 Z \beta \exp {\left( - \frac{G^*}{kT} \right)}}{N_v} (r_c - r_p)\Delta t \right)^3`


:math:`\dot{r_p} = \frac{D^r}{r_p} \frac{c^r-c^r_{eq}}{c_p^r-c^r_{eq}} + \frac{N_0 Z \beta \exp {\left( - \frac{G^*}{kT} \right)}}{N_v} (r_c - r_p)`

Now, these evolution expressions can be used to evaluate the terms :math:`L_d`, :math:`L_p` and :math:`L_s` that go straight into the CRSS expression above.




The two expressions :math:`\dot{L_p}` and :math:`\dot{L_s}` are related to each other by the common variable :math:`\dot{f_v}`. :math:`\dot{c}` is related to :math:`\dot{f_v}` as shown by the following expression



:math:`L_d = \int\dot{L_d} dt =\frac{-L_d^4}{4} \left( J_1\dot{\gamma_i} + J_2 \sum_{j \neq i} \dot{\gamma_j}  \right) - \frac{K}{2L_{d}^2}`


:math:`L_p = \left(\dot{r_p}\sqrt{\frac{2\pi}{3f_v}} - r_p\dot{f_v}\sqrt{\frac{\pi}{6f_v^3}}\right)dt`
:math:`L_s = \left(-\frac{\dot{c}}{2}\sqrt{\frac{1}{bc^3}}\right)dt`


:math:`L_d = \left(\pm \frac{\sqrt{2(KC_1 + 2})-2}{C_1}  \right)^{1/3}; C_1 = \left(J_1\dot{\gamma_i} + J_2 \sum_{j \neq i} \dot{\gamma_j}\right)`
