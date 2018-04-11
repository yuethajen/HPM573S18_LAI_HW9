import Q4_ParameterClasses as P
import Q4_MarkovModelClasses as MarkovCls
import Q4_SupportMarkovModel as SupportMarkov
import scr.SamplePathClasses as PathCls
import scr.FigureSupport as Figs

# create a cohort
cohort = MarkovCls.Cohort(
    id=1,
    therapy=P.Therapies.MONO)

# simulate the cohort
simOutputs = cohort.simulate()

# graph survival curve
PathCls.graph_sample_path(
    sample_path=simOutputs.get_survival_curve(),
    title='Survival curve',
    x_label='Simulation time step',
    y_label='Number of alive patients'
    )

# graph histogram of survival times
Figs.graph_histogram(
    data=simOutputs.get_survival_times(),
    title='Survival times of patients with atrial fibrillation',
    x_label='Survival time (years)',
    y_label='Counts',
    bin_width=1
)

# graph histogram of number of stroke
Figs.graph_histogram(
    data=simOutputs.get_survival_times(),
    title='Survival times of patients with atrial fibrillation',
    x_label='Survival time (years)',
    y_label='Counts',
    bin_width=1
)


# print the outcomes of this simulated cohort
SupportMarkov.print_outcomes(simOutputs, ' After anticoagulation:')

