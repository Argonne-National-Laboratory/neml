Yaguchi & Takahashi viscoplastic model for Grade 91 steel
=========================================================

.. warning::
   This model has been depreciated.  It produces inconistent results
   because of the time-dependent parts of the formulation.  Use at your
   own risk!

Overview
--------

This flow rule implements the complete Yaguchi & Takahashi model for
Grade 91 steel defined in [YT2000]_ and [YT2005]_.
The model has a modified Chaboche form.
This object provides a complete implementation of the flow and hardening 
functions.
Furthermore, the model hard codes the complicated interpolation formula
the original authors provide for the model coefficients.
So this implementation takes no parameters, but is nevertheless valid
over the range 473 K to 873 K.
This is the only model in NEML where units are hard-coded into the formulation,
rather than being provided by the user.

The following equations define the model:

.. math::

   \dot{\gamma}=\left\langle \frac{J_{2}\left(\operatorname{dev}\left(\bm{\sigma}\right)-\operatorname{dev}\left(\mathbf{X}\right)\right)-\sigma_{a}}{D}\right\rangle ^{n}

   J_{2}\left(\mathbf{Y}\right)=\sqrt{\frac{3}{2}\mathbf{Y}:\mathbf{Y}}

   \mathbf{g}_{\gamma}=\frac{3}{2}\frac{\operatorname{dev}\left(\bm{\sigma}\right)-\operatorname{dev}\left(\mathbf{X}\right)}{J_{2}\left(\operatorname{dev}\left(\bm{\sigma}\right)-\operatorname{dev}\left(\mathbf{X}\right)\right)}

   \bm{\alpha}=\left[\begin{array}{cccc}
   \mathbf{X}_{1} & \mathbf{X}_{2} & Q & \sigma_{a}\end{array}\right]

   \mathbf{h}_{\gamma}=\left[\begin{array}{cccc}
   \boldsymbol{X}_{1,\gamma} & \boldsymbol{X}_{2,\gamma} & Q_{\gamma} & \sigma_{a,\gamma}\end{array}\right]

   \mathbf{X}=\mathbf{X}_{1}+\mathbf{X}_{2}
   
   \mathbf{X}_{1,\gamma}=C_{1}\left(\frac{2}{3}\left(a_{10}-Q\right)\mathbf{n}-\mathbf{X}_{1}\right)\dot{\gamma}-\gamma_{1}J_{2}\left(\mathbf{X}_{1}\right)^{m-1}\mathbf{X}_{1}
   
   \mathbf{X}_{2,\gamma}=C_{2}\left(\frac{2}{3}a_{2}\mathbf{n}-\mathbf{X}_{2}\right)\dot{\gamma}-\gamma_{2}J_{2}\left(\mathbf{X}_{2}\right)^{m-1}\mathbf{X}_{2}
   
   Q_{\gamma}=d\left(q-Q\right)\sigma_{a,\gamma}	=	b\left(\sigma_{as}-\sigma_{a}\right)
   
   b	=	\begin{cases}
   b_{h} & \sigma_{as}-\sigma_{a}\ge0\\
   b_{r} & \sigma_{as}-\sigma_{a}<0
   
   \end{cases}

   \sigma_{as}	=	\left\langle A+B\log_{10}\dot{p}\right\rangle 

Parameters
----------

None, all parameters are hard coded into the object.

Class description
-----------------

.. doxygenclass:: neml::YaguchiGr91FlowRule
   :members:
   :undoc-members:
