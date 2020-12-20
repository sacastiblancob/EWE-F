# EWE-F
Modified EWE-F for Telemac coupling purposes

# Errores aparentes

- La lecutura de Spatial Ditribution estaba mal escrita, ya que la carpeta contenedora debia tener 8 caracteres para que pudiesen ser leidos los archivos de distribucion, ademas la cantidad
  de grupos funcionales estaba hardcoded a 22.
    --Arreglo: Modificar rutina readSpatialDistribution_io para que funcione en cualquier folder y con cualquier cantidad de grupos funcionales

- La variable RLEN, que define la precision de las variables de tipo real, estaba definida con un flag de compilador (#), _isDP_, se dejo hardcoded al standard de fortran para precision doble.
  Issue: Los resultados de ecosim, y por tanto de ecospace, dan completamente diferentes si se deja RLEN en precision sencilla o doble.
    --Arreglo: Dejar la precision doble, en vista de que es mejor usar variables que usen esta precision.

- En el script de ecosim.F90 se hace una llamada a la subrutina derivs.F90, el argumento "time" de la funcion se pasa harcoded en dos ocasiones como 0. y 1., cuando se esta inicializando todo,
  al activar forzado de nutrientes y produccion primaria, debido a que se esta usando precision doble, lo anterior genera un error puesto que cuando entran a la rutina se reasignan de forma extraña,
  y el 0. o 1. toman valores o muy cercanos a cero o muy grandes.
    --Arreglo: Cambiar el 0. y el 1. a 0D0 y 1D0 para dejarlos explicitamente en precision doble, y dentro de la rutina sean efectivamente o 1D0 o 0D0.

- En alguna rutina (updateForaginfTimes.F90 y ecosapce.F90) se estaba modificando la precision de una variable, a la fuerza, a precision simple real*4 (warning al momento de compilar).
    --Arreglo: se cambio para que le asignara la precision de esa variable a RLEN

- En cuanto al forzado de nutrientes y de produccion primaria se identifico que tales valores son expresados relativos a la biomasa, es decir, si las entradas de las funciones de forzado son equivalentes a 1D0
  significa que no se esta forzando nada y por tanto en ese caso la simulacion sin forzado y con forzado da igual. 1.5D por ejemplo, significaria que se esta forzando un 50% adicional de biomasa, bien a nutrientes
  o bien a los productores que se definan en los respectivos archivos de forzado de nutrientes o de forzado de produccion primaria.

# Cambios en Telemac-Coupling
- Se quitaron (comentaron) todos los condicionales para compilacion (#is_withBFM, #_isdefEcospace_, #_isdef_ForcingNutrients_, #_isdef_PrimaryProdForc_) y por tanto los procesos que se realizaban es los
  if's respectivamente, se dejo por default con ecospace activo, y para saber si se quiere usar Forzado de nutrientes o Forzado de produccion primaria se crearon dos variables en statevartypesecosim.F90
  que son tipo LOGICAL, boolFN (para Forzado de nutrientes) y boolFPP (para forzado de produccion primaria).
  Issue: Esta version solo ejecuta la version ecospace

- Makefile-ORIG tiene el makefile original de EWE-F, Makefile-ECOSPACE es el makefile para ejecutar ecospace directamente sin flags de compilacion, Makefile-MAIN es el makefie para ejecutar ecospace desde el archivo
  principal de programa main.f90 sin usar el principal ecosim.f90, tal main.f90 inicializa ecosim/ecospace y ejecuta el loop de tiempo alli mismo, el objetivo de main.f90 es poder transcribirlo directamente a
  Telemac2D.f para hacer los acoples.

- Se crea el archivo y la subrutina init_ecosim, que inicializa el programa ecosim/ecospace y prepara todo para entrar al loop temporal.

- Se crea el archivo y la subrutina time_ecosim, que realiza los calculos del loop de tiempo de ecosim

- Se crea el archivo y la subrutina clossing_and_killing, que cierra los archivos de texto abiertos y desasigna las variables de trabajo

- Muchas variables que se definian en ecospace.F90 originalmente, fueron movidas a los modulos de ecospace o ecosim, de acuerdo al caracter de tales variables.

- Se crea el archivo Makefile-MAIN-TIME, que es el makefile que se debe usar para ejectar el programa ecospace desde el programa main.f90 y no desde el programa ecosim.f90

- El programa main.f90 y sus subrutinas y modulos ya pueden usarse para acoples sincronicos con telemac2D o con cualquier otro software de modelación dinámica.

- MAKEFILE-LIB contiene la instruccion para compilar la libreria ecospace.a con "make all" que contiene los objetos necesarios para poder llamar las funciones de ecospace desde cualquier otro programa,
  en nuestro caso de interes, Telemac

!!!!!!!!!!!!!!!!! ACTUALIZACIÓN SEPTIEMBRE 18 - 2020 !!!!!!!!!!!!!!!!!!!!

- Todos los problemas aparantes de ecopath-f, ecosim-f y ecospace-f se solucionaron, después de comunicarme con el profesor Ekin, efectivamente los avaló y los dejó actualizados en los códigos fuente
  de EwE-f

- Se arregló el problema de compatibilidad entre variables nuevas y variables de WAQTEL, existía un cruce entre la posición de las variables para imprimir en VARSOR y los nombres de las variables imprimibles,
  por lo que cuando se usaba el módulo de WAQTEL, los nombres se corrían y se volvía un lío la impresión de variables en el selafin

- Corren telemac y EwE-f en simultaneo, se calcula la idoneidad de habitat, y se pueden escribir ya salidas de EwE-f en el selafin de resultados, todo lo anterior en conjunto con uso de WAQTEL

- Próximos pasos van encaminados a implementar el acople

- Hasta aquí el FortranFolder7 es funcional con todo lo anterior implementado

!!!!!!!!!!!!!!!!! ACTUALIZACIÓN DICIEMBRE 20 - 2020 !!!!!!!!!!!!!!!!!!!!

- Se decide implementar el modelo de Ecospace directamente en la malla de Telemac, para ello se añadirá como un trazador y se usará una estrategia de espacialización acorde con el compotamiento de los peces,

- El plancton definitivamente se modelará en conjunto con Waqtel (ya se incluyo el Zooplancton con el modelo del paper de la Ciénaga Grande de Santa Marta)

- El necton se modelará exclusivamente con el modelo de ecospace, si quedan acoplados estos grupos podran afectar los grupos planctónicos y los nutrientes directamente en las ecuaciones.

- Se usarán las subrutinas desarrolladas por el profesor Ekin para calcular los términos de mano derecha de las ecuaciones diferenciales.

