from dateutil.relativedelta import relativedelta

start_date = datetime.datetime(2012, 1, 1)


def stage_pars(stage):
	return


def load(data_folder, model_projection):

	riv = {}
    
    for month in range(12):
        first_of_month = start_date + relativedelta(months=month)
        
        day_of_year = first_of_month.timetuple().tm_yday
        
        stages = []
        for stage in stages:
        	lay, row, col, stage, cond, rbot = stage_pars(stage)
        	stages.append([lay, row, col, stage, cond, rbot])
        
        riv[day_of_year] = stages

	return {"stress_period_data": riv}