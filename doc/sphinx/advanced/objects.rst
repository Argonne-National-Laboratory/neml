Object and input management
===========================

Using the object system
-----------------------

NEML uses a factory system tied to a ParameterSet class to create objects
and manage input from python or XML data files.
The system is setup to incur minimal work when setting up a new object while
ensuring that the input syntax for the XML files is generated automatically
and will be consistent with the C++ library.

All a user needs to do to include a new object in the global factory and
make it available from the XML input is:

   1. Inherit from :cpp:class:`neml::NEMLObject`.
   2. Implement three static class methods:
         1. ``static std::string type()``
         2. ``static ParameterSet parameters()``
         3. ``static std::unique_ptr<NEMLObject> initialize(ParamemterSet & params)``
   3. In the header file describing the object include ``static Register<ObjectName> regObjectName`` where ``ObjectName`` is the name of the new class.

``type()`` must return a string type of the object, used to refer to it
in the factory.
This must be unique across all the NEML objects and we require it to match
the class name for classes to be merged into the main NEML project.
This example shows the implementation of this method for :cpp:class:`neml::SmallStrainPerfectPlasticity`

.. code-block:: c++

   std::string SmallStrainPerfectPlasticity::type()
   {
     return "SmallStrainPerfectPlasticity";
   }

``parameters()`` returns a ParameterSet object that informs the object system 
what parameters the object will require and sets values for default parameters,
if the object has any.  The implementation for :cpp:class:`neml::SmallStrainPerfectPlasticity` is

.. code-block:: c++

   ParameterSet SmallStrainPerfectPlasticity::parameters()
   {
     ParameterSet pset(SmallStrainPerfectPlasticity::type());

     pset.add_parameter<NEMLObject>("elastic");
     pset.add_parameter<NEMLObject>("surface");
     pset.add_parameter<NEMLObject>("ys");

     pset.add_optional_parameter<NEMLObject>("alpha",
                                             std::make_shared<ConstantInterpolate>(0.0));
     pset.add_optional_parameter<double>("tol", 1.0e-8);
     pset.add_optional_parameter<int>("miter", 50);
     pset.add_optional_parameter<bool>("verbose", false);
     pset.add_optional_parameter<int>("max_divide", 8);

     return pset;
   }

Note that other NEML objects must be added as type ``NEMLObject`` and not their
subclass type.

``initialize(ParameterSet &)`` initializes an object from a
parameter set, returning it as a ``std::unique_ptr<NEMLObject>``.
:cpp:class:`neml::SmallStrainPerfectPlasticity` implements it as

.. code-block:: c++

   std::unique_ptr<NEMLObject> SmallStrainPerfectPlasticity::initialize(ParameterSet & params)
   {
     return make_unique<SmallStrainPerfectPlasticity>(
         params.get_object_parameter<LinearElasticModel>("elastic"),
         params.get_object_parameter<YieldSurface>("surface"),
         params.get_object_parameter<Interpolate>("ys"),
         params.get_object_parameter<Interpolate>("alpha"),
         params.get_parameter<double>("tol"),
         params.get_parameter<int>("miter"),
         params.get_parameter<bool>("verbose"),
         params.get_parameter<int>("max_divide")
         ); 
   }

Finally, including the static registration class in the object header
registers it automatically with the factory.

XML input
---------

The XML input system automatically picks up the correct parameters from 
the ``parameters()`` static method.  For example, the following 
XML creates a :cpp:class:`neml::SmallStrainPerfectPlasticity` object.
Notice how the parameter names match those provided in the definition
of the C++ object.
The input system knows how to recursively ask for parameters to define
the other types of NEML objects needed.
For example, this :cpp:class:`neml::SmallStrainPerfectPlasticity` requires
a :cpp:class:`neml::LinearElasticModel` and the parameters all require
:cpp:class:`neml::Interpolate` objects.

.. code-block:: xml

   <test_perfect type="SmallStrainPerfectPlasticity">
    <elastic type="IsotropicLinearElasticModel">
      <m1 type="PolynomialInterpolate">
        <coefs>
          -100.0 100000.0
        </coefs>
      </m1>
      <m1_type>youngs</m1_type>
      <m2>0.3</m2>
      <m2_type>poissons</m2_type>
    </elastic>

    <surface type="IsoJ2"/>

    <ys type="PiecewiseLinearInterpolate">
      <points>100.0   300.0 500.0 700.0</points>
      <values>1000.0  120.0 60.0  30.0 </values>
    </ys>

    <alpha type="ConstantInterpolate">
      <v>0.1</v>
    </alpha>
   </test_perfect>


NEMLObject
----------

.. doxygenclass:: neml::NEMLObject
   :members:

ParameterSet
------------

.. doxygenclass:: neml::ParameterSet
   :members:
