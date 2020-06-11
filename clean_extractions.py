"""
    Created by: Eric March Vila (eric.march@upf.edu)
    On: 11/06/2020, 13:00 PM
"""

import pandas as pd
import pubchempy as pcp

class Cleaner(object):
    """
        Class with functions for cleaning eTox data extractions from Vitic database
    """

    def __init__(self, dataframe: pd.DataFrame, timepoint_field: str,
                                                timepoint_unit_field: str,
                                                animal_age_field: str,
                                                animal_age_unit_field: str,
                                                parameter_field: str,
                                                value_field: float,
                                                value_unit_field: str):
                                            
        self.df = dataframe
        self.timepoint_field = timepoint_field
        self.timepoint_unit_field = timepoint_unit_field
        self.animal_age_field = animal_age_field
        self.animal_age_unit_field = animal_age_unit_field
        self.parameter_field = parameter_field
        self.value_field = value_field
        self.value_unit_field = value_unit_field

    def clean_dataframe(self) -> pd.DataFrame:
        """
            Main function. It runs the other functions to clean the dataframe:
                - removing NaN values
                - converting units to international unit system. Ex: mg/100mL or mg% to mM
                - converting all time units from weeks, months and years to days

            :return cleaned_dataframe:
        """
        
        # Remove NaN and negative age and timepoint
        self.remove_nans()

        # Changes animal age units into days
        self.animal_age_transformation()

        # Changes timepoint units into days
        self.timepoint_transformation()

        # Corrects the units. Uses the weight of the parameter (i.e. Glucose) to get the proper concentration
        
        concentration_units = self.get_concentration_units()
        param_from_conc_unit_list = self.df.loc[self.df[self.value_unit_field].isin(concentration_units), self.parameter_field].unique()
        parameter_weight = self.get_parameter_weights(param_from_conc_unit_list)
        cleaned_dataframe = self.correct_units(concentration_units, parameter_weight)

        return cleaned_dataframe

    def remove_nans(self):
        """
            Removes NaN and negative values from Timepoint and age
        """

        self.df.drop(self.df[self.df[self.timepoint_field] < 0].index, axis=0, inplace=True)
        self.df.dropna(subset=[self.timepoint_field,self.animal_age_field], inplace=True)
    
    def animal_age_transformation(self):
        """
            Transforms animal age into days
        """

        age_start_days_from_weeks = self.df.loc[self.df[self.animal_age_unit_field].isin(['Weeks', 'Week', 'weeks', 'Weels', 'Weeka']), self.animal_age_field] * 7
        age_start_days_from_months = self.df.loc[self.df[self.animal_age_unit_field].isin(['Months','Month','Moths','Monhts']), self.animal_age_field] * 30
        age_start_days_from_years = self.df.loc[self.df[self.animal_age_unit_field].isin(['Years','Year']), self.animal_age_field] * 360
        age_start_days_from_other_units = self.df.loc[self.df[self.animal_age_unit_field].isin(['Days', 'Day','0',]), self.animal_age_field].values

        self.df.loc[self.df[self.animal_age_unit_field].isin(['Weeks', 'Week', 'weeks', 'Weels', 'Weeka']), 'age_at_start(days)'] = age_start_days_from_weeks
        self.df.loc[self.df[self.animal_age_unit_field].isin(['Months','Month','Moths','Monhts']), 'age_at_start(days)'] = age_start_days_from_months
        self.df.loc[self.df[self.animal_age_unit_field].isin(['Years','Year']), 'age_at_start(days)'] = age_start_days_from_years
        self.df.loc[self.df[self.animal_age_unit_field].isin(['Days', 'Day','0',]), 'age_at_start(days)'] = age_start_days_from_other_units
    
    def timepoint_transformation(self):
        """
            Transforms timepoint units into days
        """

        timepoint_days_from_weeks = self.df.loc[self.df[self.timepoint_unit_field].isin(['Weeks', 'Week','WEEK','week']), self.timepoint_field] * 7
        timepoint_days_from_months = self.df.loc[self.df[self.timepoint_unit_field].isin(['Months','Month']), self.timepoint_field] * 30
        timepoint_days_from_years = self.df.loc[self.df[self.timepoint_unit_field].isin(['Years']), self.timepoint_field] * 360
        timepoint_days_from_hours = self.df.loc[self.df[self.timepoint_unit_field].isin(['Hours']), self.timepoint_field] / 24
        timepoint_days_from_minutes = self.df.loc[self.df[self.timepoint_unit_field].isin(['Minutes']), self.timepoint_field] / 1440
        timepoint_days_from_other_units = self.df.loc[self.df[self.timepoint_unit_field].isin(['Days','Day','Day(s)','day']), self.timepoint_field].values

        self.df.loc[self.df[self.timepoint_unit_field].isin(['Weeks', 'Week','WEEK','week']), 'timepoint(days)'] = timepoint_days_from_weeks
        self.df.loc[self.df[self.timepoint_unit_field].isin(['Months','Month']), 'timepoint(days)'] = timepoint_days_from_months
        self.df.loc[self.df[self.timepoint_unit_field].isin(['Years']), 'timepoint(days)'] = timepoint_days_from_years
        self.df.loc[self.df[self.timepoint_unit_field].isin(['Hours']), 'timepoint(days)'] = timepoint_days_from_hours
        self.df.loc[self.df[self.timepoint_unit_field].isin(['Minutes']), 'timepoint(days)'] = timepoint_days_from_minutes
        self.df.loc[self.df[self.timepoint_unit_field].isin(['Days','Day','Day(s)','day']), 'timepoint(days)'] = timepoint_days_from_other_units

        ### Sums the timepoint and the age to get the age at the end of the experiment
        self.df['timepoint_age(days)'] = self.df.iloc[:,-2:].sum(axis=1)
    
    def get_concentration_units(self) -> list:
        """
            Gets the units for the average value.

            :return conc_units:
        """

        conc_units = []
        for unit in self.df[self.value_unit_field].unique():
            if 'mg' in unit.lower():
                conc_units.append(unit)
            elif unit.lower().strip('(').startswith('g'):
                conc_units.append(unit)
            elif 'mol' in unit.lower():
                conc_units.append(unit)
            elif 'milligrams' in unit.lower():
                conc_units.append(unit)
            elif unit.lower().startswith('ng'):
                conc_units.append(unit)
            elif unit.lower().startswith('ug'):
                conc_units.append(unit)
            elif 'pg' in unit.lower():
                conc_units.append(unit)
            elif 'ratio' in unit.lower():
                conc_units.append(unit)
            elif unit.lower().startswith('mm'):
                conc_units.append(unit)
            elif unit.lower().startswith('mcg'):
                conc_units.append(unit)
            elif 'm g ' in unit.lower():
                conc_units.append(unit)
            elif unit.lower().startswith('mcm'):
                conc_units.append(unit)
            elif 'microg' in unit.lower():
                conc_units.append(unit)
            elif unit.lower().startswith('um'):
                conc_units.append(unit)
            elif unit.lower().startswith('nm'):
                conc_units.append(unit)
            elif 'pmoi' in unit.lower():
                conc_units.append(unit)
            elif unit.startswith('MNOL'):
                conc_units.append(unit)
            elif 'nrnol/L' in unit:
                conc_units.append(unit)
            elif 'm e q' in unit.lower() or 'meq' in unit.lower():
                conc_units.append(unit)
            elif 'MAEQ/L' in unit or 'ME Q/L' in unit:
                conc_units.append(unit)
        
        return conc_units
    
    def get_parameter_weights(self, parameter_conc_unit_list: list) -> dict:
        """
            Get the molecular weight of the selected parameter.

            :param parameter_conc_unit_list: list of concentration units of the parameter

            :return parameter_weight:
        """

        parameter_weight = {}
        
        for parameter in parameter_conc_unit_list:
            compounds = pcp.get_compounds(parameter, 'name')
            if compounds:
                for compound in compounds:
                    mol_weight = compound.molecular_weight
                    mol_formula = compound.molecular_formula
                    if mol_weight:
                        parameter_weight.update({parameter:mol_weight})
        
        return parameter_weight
            
    def convert_values_to_SI(value: float, weight: float, arguments: list) -> float: 
        """
            Function to convert the values to the International Unit System

            :param value: old value
            :param weight: weight of the molecule in order to obtain its molar concentration
            :arguments: list of units to be transformed

            :return standard_val: standardized value into IUS
        """

        if 'micro' in arguments:
            standard_val = value*1e-3
        elif 'nano' in arguments:
            standard_val = value*1e-6
        elif 'pico' in arguments:
            standard_val = value*1e-9
        else:
            standard_val = value
        
        if 'gram' in arguments:
            standard_val = (standard_val*1000)/weight
        elif 'mol' in arguments:
            standard_val = standard_val*1000
        elif 'miligram' in arguments:
            standard_val = standard_val/weight
            
        if 'mililiter' in arguments:
            standard_val = standard_val*1e-3
        elif 'deciliter' in arguments:
            standard_val = standard_val*1e-1
        
        if 'minute' in arguments:
            standard_val = standard_val*60
            
        return standard_val

    def values_to_SI(param_per_unit_df: pd.DataFrame, parameter_weight: dict, *args) -> list:
        """
            Checks each value amd converts it into its IUS equivalent

            :param param_per_unit_df: dataframe containing the parameter and the associated average value
            :param *args: list of units to be converted

            :return list_of_vals: list of converted values
        """

        list_of_vals = []
        print(param_per_unit_df)
        for index, row in param_per_unit_df.iterrows():
            param = row[0]
            value = row[1]
            if param in parameter_weight.keys():
                weight = parameter_weight[param]
                std_val = convert_values_to_SI(value,weight,args[0])
                list_of_vals.append(std_val)
                
        return list_of_vals
    
    def correct_units(self, concentration_units: list, parameter_weight: dict) -> pd.DataFrame:
        """
            Function that applies the correction to the units.

            :param concentration_units:
            :param parameter_weight:

            :return corrected_df:
        """

        for unit in concentration_units:
            unit_split = unit.split('/')
            args_ = []

            if  len(unit_split) == 2:
                unit_, liter = [u.lower() for u in unit_split]
                
                if 'u' in unit_ or unit_.startswith('mc') or 'micro' in unit_:
                    args_.append('micro')
                elif unit_.startswith('p'):
                    args_.append('pico')
                elif unit_.startswith('n'):
                    args_.append('nano')
                
                if unit_.lower().strip('?').strip('*').strip('<').startswith('mol'):
                    args_.append('mol')
                if unit_.replace(' ','').replace('gm','mg').startswith('mg') or unit_.startswith('mk') or 'milligrams' in unit_.lower():
                    args_.append('miligram')
                elif unit_.lower().strip('(').startswith('g'):
                    args_.append('gram')
                
                if 'dl' in liter or 'deciliter' in liter or '100ml' in liter.replace(' ','').lower():
                    args_.append('deciliter')
                elif 'ml' in liter:
                    args_.append('mililiter')
                    
            elif len(unit_split) > 2:
                unit_, time, liter = [u.lower() for u in unit_split]
                if 'n' in unit_:
                    args_.append('nano')
                elif 'u' in unit_:
                    args_.append('micro')
                
                if 'min' in time:
                    args_.append('minute')
                elif 'ml' in time:
                    args_.append('mililiter')
                    
                if 'min' in liter:
                    args_.append('minute')
                elif liter.startswith('g'):
                    args_.append('gram')

            else:
                unit_ = unit_split[0].lower()
                if unit_.startswith('n'):
                    args_.append('nano')
                elif unit_.startswith('u') or 'micro' in unit_:
                    args_.append('micro')
                    
                if unit_.startswith('g'):
                    args_.append('gram')
                elif unit_.startswith('mg'):
                    args_.append('miligram')
                    
                if '%' in unit_ or '100ml' in unit_.replace(' ','').lower():
                    args_.append('deciliter')
            
            param_per_unit = self.df.loc[(self.df[self.value_unit_field] == unit),
                                                    [self.parameter_field,self.value_field]]
            
            list_of_vals = self.values_to_SI(param_per_unit, parameter_weight, args_)
            
            self.df.loc[(self.df[self.value_unit_field] == unit)  & (self.df[self.parameter_field].isin(parameter_weight.keys())),
                                    'average_value_fixed'] = list_of_vals
            self.df.loc[(self.df[self.value_unit_field] == unit),'unit_fixed'] = 'mM'
        
        corrected_df = self.df.copy()

        return corrected_df