# Config Variables
CONFIG = dict (
		InitialRecessiveAlleleFraction = 0.5,     # the fraction of short-lived population at start
		equilibrium_population_size = dict(),
		max_experimental_age = dict(),
		avg_lifespan = dict(),
		age_distribution = dict(		  # This is for storing age distribution data for your current type of run
			CDF = dict(),
     		Absolute = dict(),
     		PDF = dict(),
       		change_per_x_unit = dict(),
			hazard = dict(),
		),								  # Should be stored in the format CONFIG['age_distribution'][<distribution type>][<population name>]
		equilibrium_lifespan_distribution = dict(		  # This is for storing equilibrium_lifespan_distribution data for your current type of run
			CDF = dict(),
     		Absolute = dict(),
     		PDF = dict(),
			change_per_x_unit = dict(),
			hazard = dict(),
    	),								  # Should be stored in the format CONFIG['equilibrium_lifespan_distribution'][<distribution type>][<population name>]
        lifespan_distribution = dict(		  # This is for storing lifespan distribution data for your current type of run
			CDF = dict(),
     		Absolute = dict(),
     		PDF = dict(),
			change_per_x_unit = dict(),
			hazard = dict(),
    	),	
          ComparisonDistibution = dict(	  # This is for when you want to graph a second, pre-run distribution with the opposite value for UseLinearGrowth
			lifespan_distribution = dict(
                CDF = dict(),
				Absolute = dict(),
				PDF = dict(),
				change_per_x_unit = dict(),
				hazard = dict(),
			),
   			equilibrium_lifespan_distribution = dict(
                CDF = dict(),
				Absolute = dict(),
				PDF = dict(),
				change_per_x_unit = dict(),
				hazard = dict(),
			),
			age_distribution = dict(
                CDF = dict(),
				Absolute = dict(),
				PDF = dict(),
			)
		)	  							  # Should be stored in the format CONFIG['ComparisonDistibution'][<data_name>][<distribution type>][<population name>]

  )

