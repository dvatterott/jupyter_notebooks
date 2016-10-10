
#basic analyses script

def dataTableRT(ValueDictionary):
    import numpy
    sortedD = sorted(ValueDictionary)
    RedString = '\x1b[1;31m'
    BlackString = '\x1b[0m'
    headerlist=[]
    print (100*'-')
    
    dataMat = numpy.zeros(shape=(len(ValueDictionary[sortedD[0]]),len(ValueDictionary.keys()))) #make a matrix to put the data into
    for headerCounter, header in enumerate(sortedD): 
        headerlist.append(header)
        for rowCounter, row in enumerate(ValueDictionary[header]): dataMat[rowCounter,headerCounter] = row
    
    stringformatHEADER = '%s\t\t\t'*(len(headerlist)+1)
    StringVariables = ['Sub#']
    for header in headerlist: StringVariables.append(header)
    StringVariables = tuple(StringVariables)
    print RedString+stringformatHEADER % StringVariables+BlackString
    
    stringformatData = '%s\t\t\t'*(len(headerlist)+1)
    for row in xrange(0,len(dataMat)):
        StringVariablesData = [str(row+1)]
        for SubMeans in dataMat[row]: StringVariablesData.append(round(SubMeans,2))
        StringVariablesData = tuple(StringVariablesData)
        print stringformatData % StringVariablesData
        
    means = numpy.mean(dataMat,axis=0)
    StringVariables = ['mean']
    for CondMean in means: StringVariables.append(round(CondMean,2))
    StringVariables = tuple(StringVariables)
    print RedString + stringformatData % StringVariables
    
    STEs = (numpy.std(dataMat,axis=0, ddof=1))/numpy.sqrt(len(dataMat)) #axis0 is collapse across rows. ddof 1 makes the std into a sample std
    StringVariables = ['STE']
    for CondSTE in STEs: StringVariables.append(round(CondSTE,2))
    StringVariables = tuple(StringVariables)
    print stringformatData % StringVariables + BlackString
    print(100*'-')

def dataTableAcc(ValueDictionary):
    import numpy
    sortedD = sorted(ValueDictionary)
    RedString = '\x1b[1;31m'
    BlackString = '\x1b[0m'
    headerlist=[]
    print (100*'-')
    
    dataMat = numpy.zeros(shape=(len(ValueDictionary[sortedD[0]]),len(ValueDictionary.keys()))) #make a matrix to put the data into
    for headerCounter, header in enumerate(sortedD): 
        headerlist.append(header)
        for rowCounter, row in enumerate(ValueDictionary[header]): dataMat[rowCounter,headerCounter] = row
    
    stringformatHEADER = '%s\t\t\t'*(len(headerlist)+1)
    StringVariables = ['Sub#']
    for header in headerlist: StringVariables.append(header)
    StringVariables = tuple(StringVariables)
    print RedString+stringformatHEADER % StringVariables+BlackString
    
    stringformatData = '%s\t\t\t'*(len(headerlist)+1)
    for row in xrange(0,len(dataMat)):
        StringVariablesData = [str(row+1)]
        for SubMeans in dataMat[row]: StringVariablesData.append(round(SubMeans,2))
        StringVariablesData = tuple(StringVariablesData)
        print stringformatData % StringVariablesData
        
    means = numpy.mean(dataMat,axis=0)
    StringVariables = ['mean']
    for CondMean in means: StringVariables.append(round(CondMean,2))
    StringVariables = tuple(StringVariables)
    print RedString + stringformatData % StringVariables
    
    ErrorRates = (1-means)*100
    StringVariables = ['Error%']
    for CondER in ErrorRates: StringVariables.append(round(CondER,2))
    StringVariables = tuple(StringVariables)
    print stringformatData % StringVariables + BlackString
    print(100*'-')

def OneWayRepeatedMeasures_ANOVA(ValueDictionary):
    import scipy.stats, numpy
    print '1 way repeated-measures ANOVA:'
    
    sortedD = sorted(ValueDictionary)
    headerlist=[]
    dataMat = numpy.zeros(shape=(len(ValueDictionary[sortedD[0]]),len(ValueDictionary.keys()))) #make a matrix to put the data into
    for headerCounter, header in enumerate(sortedD): 
        headerlist.append(header)
        for rowCounter, row in enumerate(ValueDictionary[header]): dataMat[rowCounter,headerCounter] = row
            
    SSmodel=0
    df_model = len(headerlist)-1
    means = numpy.mean(dataMat,axis=0)
    for mean in means: SSmodel+=len(dataMat)*numpy.square(mean-numpy.mean(dataMat))
    SSwithin=0
    df_within = len(dataMat)*df_model
    for rows in dataMat: SSwithin+=(numpy.var(rows, ddof=1)*(len(rows)-1))
    SSresidual = SSwithin-SSmodel
    df_residual = df_within - df_model
    MeanSquare_model = SSmodel/df_model
    MeanSquare_residual = SSresidual/df_residual
    F_value = MeanSquare_model/MeanSquare_residual
    p_value = scipy.stats.f.sf(F_value,df_model,df_within) #note this should actually be corrected for violation of sphericity
    print '\tF value = %s' % (str(round(F_value,2))) + '\tp value = %s' % (str(round(p_value,4)))

def OneWayConfInterval(table):
    import scipy.stats, numpy
    ParticipantsMeans = []
    STEs = []
    CIs = []
    for participant in xrange(numpy.shape(table)[0]):
        mean = []
        for condition in xrange(numpy.shape(table)[1]):
            mean.append(table[condition][participant+1])
        ParticipantsMeans.append(sum(mean)/len(mean))
    ConfMeans = numpy.zeros(shape=numpy.shape(table))
    for participant in xrange(numpy.shape(table)[0]):
        for condition in xrange(numpy.shape(table)[1]):
            ConfMeans[participant][condition] = table[condition][participant+1]-\
            ParticipantsMeans[participant]+numpy.array(ParticipantsMeans).mean()
    for counter, column in enumerate(ConfMeans.T):
        STEs.append(numpy.std(column, ddof=1)/numpy.sqrt(len(column)))
        CIs.append(STEs[counter]*scipy.stats.t.isf(0.025, len(ConfMeans)-1))
    return CIs
    
def SimpleComparisonCI(table):
    import scipy.stats, math
    ttest = scipy.stats.ttest_rel(table[0], table[1])
    MeanDiff_byT = abs((table[0].mean()-table[1].mean())/ttest[0])
    CI = MeanDiff_byT*scipy.stats.t.isf(0.025, len(table)-1)*math.pow(2,0.05)/2
    return CI

def SimpleComparisonCI_1(table):
    import scipy.stats, math
    ttest = scipy.stats.ttest_rel(table[1], table[2])
    MeanDiff_byT = abs((table[1].mean()-table[2].mean())/ttest[0])
    CI = MeanDiff_byT*scipy.stats.t.isf(0.025, len(table)-1)*math.pow(2,0.05)/2
    return CI

def Fit_PowerFunction2(x, y):
    ##########
    # Fitting the data -- Least Squares Method
    ##########

    # Power-law fitting is best done by first converting
    # to a linear equation and then fitting to a straight line.
    #
    #  y = a * x^b
    #  log(y) = log(a) + b*log(x)
    # this is from wiki.scipy.org/Cookbook/FittingData
    import scipy, numpy
    logx = scipy.log10(x); logy = scipy.log10(y) #log the data so that we can use a linear func
    fitfunc = lambda p, x: p[0] + p[1] * x #target fun. notice this is a linear fun. 
    errfunc = lambda p, x, y: (y - fitfunc(p, x)) #distance from our target fun
    pinit = [750.0, -0.05] #inital parameter guess
    pfinal,covar,infodict,mesg,ier = scipy.optimize.leastsq(errfunc, pinit, args=(logx, logy), full_output=1)    
    index = pfinal[1]
    amp = 10.0**pfinal[0]
    ss_err = (infodict['fvec']**2).sum()
    ss_tot = ((y-y.mean())**2).sum()
    rsquared = 1-(ss_err/ss_tot)
    return (pfinal, covar, index, amp, rsquared)

def Fit_ExpDecFunction2(x, y):
    ##########
    # Fitting the data -- Least Squares Method
    ##########
    # Exp Decay fitting
    import scipy, numpy as np
    fitfunc = lambda p, x: p[0]*np.exp(-p[1]*x) + p[2]#target fun. notice this is a linear fun. 
    errfunc = lambda p, x, y: (y - fitfunc(p, x)) #distance from our target fun
    pinit = [50.0, 0.1, 850] #inital parameter guess
    pfinal,covar,infodict,mesg,ier = scipy.optimize.leastsq(errfunc, pinit, args=(x, y), full_output=1)    
    index = pfinal[1]
    amp = pfinal[0]
    asmp = pfinal[2]
    ss_err = (infodict['fvec']**2).sum()
    ss_tot = ((y-y.mean())**2).sum()
    rsquared = 1-(ss_err/ss_tot)
    dof = len(x)-len(pfinal)
    rmse = np.sqrt(ss_err/dof)
    return (rmse, index, amp, asmp, rsquared)

def Fit_ExpDecFunction(x, y, bounders, guess):
    ##########
    # Fitting the data -- Least Squares Method
    ##########
    # Exp Decay fitting
    import scipy, numpy as np, leastsqbound
    fitfunc = lambda p, x: p[0]*np.exp(-p[1]*x) + p[2]#target fun. notice this is a linear fun. 
    errfunc = lambda p, x, y: (y - fitfunc(p, x)) #distance from our target fun
    pinit = guess #inital parameter guess
    pfinal,covar,infodict,mesg,ier = leastsqbound.leastsqbound(errfunc, pinit, args=(x, y), bounds=bounders,full_output=1)    
    index = pfinal[1]
    amp = pfinal[0]
    asmp = pfinal[2]
    ss_err = (infodict['fvec']**2).sum()
    ss_tot = ((y-y.mean())**2).sum()
    rsquared = 1-(ss_err/ss_tot)
    dof = len(x)-len(pfinal)
    rmse = np.sqrt(ss_err/dof)
    return (rmse, index, amp, asmp, rsquared)

def Fit_PowerFunction_OldEq(x,y):
    # Power-law fitting is best done by first converting
    # to a linear equation and then fitting to a straight line.
    #
    #  y = a * x^b
    #  log(y) = log(a) + b*log(x)
    # this is from wiki.scipy.org/Cookbook/FittingData, but with lin reg
    import scipy, numpy
    logx = scipy.log10(x); logy = scipy.log10(y) #log the data so that we can use a linear func
    slope, intercept, r_value, p_value, stderr = scipy.stats.linregress(logx,logy)
    amp = 10.0**intercept
    index = slope
    rsquared = r_value**2
    return (amp,index,rsquared,p_value,stderr)

def Fit_PowerFunction(x,y):
    # Power-law fitting is best done by first converting
    # to a linear equation and then fitting to a straight line.
    #
    #  y = a * x^b
    #  log(y) = log(a) + b*log(x)
    # this is from wiki.scipy.org/Cookbook/FittingData, but with lin reg
    import scipy, numpy as np, math
    fitfunc = lambda p, x: p[0]*(x**p[1]) + p[2]#target fun. notice this is a linear fun. 
    errfunc = lambda p, x, y: (y - fitfunc(p, x)) #distance from our target fun
    pinit = (1000, -0.5, 800) #inital parameter guess
    pfinal,covar,infodict,mesg,ier = scipy.optimize.leastsq(errfunc, pinit, args=(x, y), full_output=1) 
    index = pfinal[1]
    amp = pfinal[0]
    asmp = pfinal[2]
    ss_err = (infodict['fvec']**2).sum()
    ss_tot = ((y-y.mean())**2).sum()
    rsquared = 1-(ss_err/ss_tot)
    dof = len(x)-len(pfinal)
    rmse = np.sqrt(ss_err/dof)
    return (rmse, index, amp, asmp, rsquared)
    
def Print_PowerFunFitOut2(condition, datain):
    print '----------------------------------------'
    print condition
    print 'index=' + str(datain[2])
    print 'amp=' + str(datain[3])
    print 'r^2=' + str(datain[4])

def Print_PowerFunFitOut(condition, datain):
    print '----------------------------------------'
    print condition
    print 'index=' + str(datain[1])
    print 'amp=' + str(datain[0])
    print 'r^2=' + str(datain[2])
    print 'p_value=' + str(datain[3])
    print 'stderr=' + str(datain[4])

def Plot_PowerFunFitOut(amp_Absent, index_Absent, amp_Present, index_Present, xdata, ydata_Absent, ydata_Present, ylim=[500, 1200]):
    import matplotlib
    import matplotlib.pyplot as plt
    yFunData_Absent = amp_Absent*(xdata**index_Absent)
    yFunData_Present = amp_Present*(xdata**index_Present)

    fig = plt.figure(figsize=(15,4))
    axes = fig.add_subplot(111)

    #normal plot
    axes.plot(xdata, ydata_Absent, 'bo', markersize=5, label="_")
    axes.plot(xdata, ydata_Present, 'gs', markersize=5, label="_")
    axes.plot(xdata,  yFunData_Absent, 'blue', label="Singleton Absent")
    axes.plot(xdata, yFunData_Present, 'green', label="Singleton Present")
    axes.legend()
    axes.set_title('Best Fit Power Law')
    axes.set_ylabel('Response Times (ms)')
    axes.set_xlabel('# Experiences')
    axes.set_xlim(1, max(xdata))
    axes.set_ylim(ylim[0], ylim[1])
    #text(5, 6.5, 'Ampli = %5.2f +/- %5.2f' % (amp, ampErr))
    #text(5, 5.5, 'Index = %5.2f +/- %5.2f' % (index, indexErr))

    #log plot
    #axes[1].loglog(xdata, ydata_Absent, 'bo', markersize=5, label="_")
    #axes[1].loglog(xdata, ydata_Present, 'gs', markersize=5, label="_")
    #axes[1].plot(xdata,  yFunData_Absent, 'blue',label="Singleton Absent")
    #axes[1].plot(xdata, yFunData_Present, 'green',label="Singleton Present")
    #axes[1].set_title('Best Fit Power Law (log scale)')
    #axes[1].set_ylabel('Response Times (log scale)')
    #axes[1].set_xlabel('# Experiences (log scale)')
    #axes[1].set_xlim(1, max(xdata))
    #axes[1].set_ylim(ylim[0], ylim[1])

def Test_and_Plot_2x2_LearnDistRej(tableRT, low, high):
    #t-tests for RT data
    import numpy as np, scipy.stats, pandas as pd, PythonAnalyses, matplotlib, matplotlib.pyplot as plt, pylab as pl
    ttest1 = scipy.stats.ttest_rel(tableRT[:][1][1], tableRT[:][1][2]) #ttests. the_rel is within subjects
    ttest2 = scipy.stats.ttest_rel(tableRT[:][2][1], tableRT[:][2][2])

    normTest1 = scipy.stats.shapiro(tableRT[:][1][1]-tableRT[:][1][2]) #shapiro-wilk test for normality
    normTest2 = scipy.stats.shapiro(tableRT[:][2][1]-tableRT[:][2][2])
    print '--------------------------------------------------------'
    print 'Normality Tests'
    print 'Test 1 ='+str(normTest1)
    print 'Test 2 ='+str(normTest2)
    print '--------------------------------------------------------'
    print 'T-Test: 1st Block Half - SingletonAbsent vs. SingletonPresent'
    print '\tt value = %s' % (str(round(abs(ttest1[0]),2))) + '\tp value = %s' % (str(round(ttest1[1],4)))
    print 'T-Test: 2nd Block Half - SingletonAbsent vs. SingletonPresent'
    print '\tt value = %s' % (str(round(abs(ttest2[0]),2))) + '\tp value = %s' % (str(round(ttest2[1],4)))

    CIs = [0,0]
    CIs[0] = PythonAnalyses.SimpleComparisonCI_1(tableRT[1])
    CIs[1] = PythonAnalyses.SimpleComparisonCI_1(tableRT[2])
    PlotFrame = pd.DataFrame([tableRT[(1)].mean(), tableRT[(2)].mean()])
    PlotFrame.columns = ['Singleton Absent', 'Singleton Present']
    PlotFrame.plot(ylim = [low, high], kind='bar', yerr=CIs)
    plt.xticks(range(2), ['1st Half of Block','2nd Half of Block'], rotation=0)
    print 'n = %s' % str(np.size(tableRT.index))

def Test_and_Plot_PairwiseComparison(df, low, high, xTicks, title):
    import numpy as np, scipy.stats, pandas as pd, PythonAnalyses, matplotlib, matplotlib.pyplot as plt, pylab as pl
    ttest1 = scipy.stats.ttest_rel(df[1], df[2]) #ttests. the_rel is within subjects
    print 'ttest = %s' % str(ttest1)

    CIs = [0]
    CIs = PythonAnalyses.SimpleComparisonCI_1(df)
    PlotFrame = pd.DataFrame([df[(1)].mean(), df[(2)].mean()])
    PlotFrame.plot(ylim = [low, high], kind='bar', yerr=CIs)
    plt.xticks(range(2), xTicks, rotation=0)
    plt.legend().set_visible(False)
    plt.title(title)
    print 'n = %s' % str(np.size(df.index))

def Test_and_Plot_PairwiseComparison2(df, low, high, xTicks, title):
    import numpy as np, scipy.stats, pandas as pd, PythonAnalyses, matplotlib, matplotlib.pyplot as plt, pylab as pl
    ttest1 = scipy.stats.ttest_rel(df[0], df[1]) #ttests. the_rel is within subjects
    print 'ttest = %s' % str(ttest1)

    CIs = [0]
    CIs = PythonAnalyses.SimpleComparisonCI(df)
    PlotFrame = pd.DataFrame([df[(0)].mean(), df[(1)].mean()])
    PlotFrame.plot(ylim = [low, high], kind='bar', yerr=CIs)
    plt.xticks(range(2), xTicks, rotation=0)
    plt.legend().set_visible(False)
    plt.title(title)
    print 'n = %s' % str(np.size(df.index))

def average_acrossBlock_PowerFunTests(tableRT):
    import pandas as pd, PythonAnalyses, numpy as np, scipy.stats
    Ab_Frame = pd.DataFrame(np.zeros([np.size(tableRT.index),5]),index=tableRT.index,columns=['amp','index','rsquared','p_value','stderr'])
    Pre_Frame = pd.DataFrame(np.zeros([np.size(tableRT.index),5]),index=tableRT.index,columns=['amp','index','rsquared','p_value','stderr'])

    for subs in tableRT.index:
        FiniteList1 = np.isfinite(tableRT.loc[subs][:,1]); FiniteList2 = np.isfinite(tableRT.loc[subs][:,2])
        AbsentOut = PythonAnalyses.Fit_PowerFunction(tableRT.loc[subs][:,1][FiniteList1].index, tableRT.loc[subs][:,1][FiniteList1])
        PresentOut = PythonAnalyses.Fit_PowerFunction(tableRT.loc[subs][:,2][FiniteList2].index, tableRT.loc[subs][:,2][FiniteList2])
        Ab_Frame.loc[subs] = pd.Series({'amp':AbsentOut[0],'index':AbsentOut[1],'rsquared':AbsentOut[2],'p_value':AbsentOut[3],'stderr':AbsentOut[4]})
        Pre_Frame.loc[subs] = pd.Series({'amp':PresentOut[0],'index':PresentOut[1],'rsquared':PresentOut[2],'p_value':PresentOut[3],'stderr':PresentOut[4]})

    ttest1 = scipy.stats.ttest_rel(Pre_Frame['index'], Ab_Frame['index'])
    ttest2 = scipy.stats.ttest_rel(Pre_Frame['amp'], Ab_Frame['amp'])

    normTest1 = scipy.stats.shapiro(Pre_Frame['index']-Ab_Frame['index']) #shapiro-wilk test for normality
    normTest3 = scipy.stats.shapiro(Pre_Frame['amp']-Ab_Frame['amp'])

    print '--------------------------------------------------------'
    print 'Normality Tests'
    print 'Test1 ='+str(normTest1)
    print 'Test2 ='+str(normTest3)
    print '--------------------------------------------------------'
    print 'index ttest='+str(ttest1)
    print 'amp ttest='+str(ttest2)
    return (Ab_Frame, Pre_Frame)

def average_acrossBlock_ExpDecFunTests(tableRT):
    import pandas as pd, PythonAnalyses, numpy as np, scipy.stats
    Ab_Frame = pd.DataFrame(np.zeros([np.size(tableRT.index),5]),index=tableRT.index,columns=['amp','index','asmp','rsquared','rmse'])
    Pre_Frame = pd.DataFrame(np.zeros([np.size(tableRT.index),5]),index=tableRT.index,columns=['amp','index','asmp','rsquared','rmse'])

    for subs in tableRT.index:
        FiniteList1 = np.isfinite(tableRT.loc[subs][:,1]); FiniteList2 = np.isfinite(tableRT.loc[subs][:,2])
        AbsentOut = PythonAnalyses.Fit_ExpDecFunction(tableRT.loc[subs][:,1][FiniteList1].index, tableRT.loc[subs][:,1][FiniteList1])
        PresentOut = PythonAnalyses.Fit_ExpDecFunction(tableRT.loc[subs][:,2][FiniteList2].index, tableRT.loc[subs][:,2][FiniteList2])
        Ab_Frame.loc[subs] = pd.Series({'amp':AbsentOut[2],'index':AbsentOut[1],'rsquared':AbsentOut[4],'asmp':AbsentOut[3],'rmse':AbsentOut[0]})
        Pre_Frame.loc[subs] = pd.Series({'amp':PresentOut[2],'index':PresentOut[1],'rsquared':PresentOut[4],'asmp':PresentOut[3],'rmse':PresentOut[0]})

    ttest1 = scipy.stats.ttest_rel(Pre_Frame['index'], Ab_Frame['index'])
    ttest2 = scipy.stats.ttest_rel(Pre_Frame['amp'], Ab_Frame['amp'])
    ttest3 = scipy.stats.ttest_rel(Pre_Frame['asmp'], Ab_Frame['asmp'])

    normTest1 = scipy.stats.shapiro(Pre_Frame['index']-Ab_Frame['index']) #shapiro-wilk test for normality
    normTest2 = scipy.stats.shapiro(Pre_Frame['amp']-Ab_Frame['amp'])
    normTest3 = scipy.stats.shapiro(Pre_Frame['asmp']-Ab_Frame['asmp'])

    print '--------------------------------------------------------'
    print 'Normality Tests'
    print 'Test1 ='+str(normTest1)
    print 'Test2 ='+str(normTest2)
    print 'Test3 ='+str(normTest3)
    print '--------------------------------------------------------'
    print 'index ttest='+str(ttest1)
    print 'amp ttest='+str(ttest2)
    print 'asmp ttest='+str(ttest3)
    return (Ab_Frame, Pre_Frame)

def Plot_ExpDecFunFitOut(AbsentFit, PresentFit, xdata, ydata_Absent, ydata_Present, ylim=[500, 1200]):
    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np
    yFunData_Absent = AbsentFit[2]*np.exp(-AbsentFit[1]*xdata)+AbsentFit[3]
    yFunData_Present = PresentFit[2]*np.exp(-PresentFit[1]*xdata)+PresentFit[3]

    fig = plt.figure(figsize=(15,4))
    axes = fig.add_subplot(111)

    #normal plot
    axes.plot(xdata, ydata_Absent, 'bo', markersize=5, label="_")
    axes.plot(xdata, ydata_Present, 'gs', markersize=5, label="_")
    axes.plot(xdata,  yFunData_Absent, 'blue', label="Singleton Absent")
    axes.plot(xdata, yFunData_Present, 'green', label="Singleton Present")
    axes.legend()
    axes.set_title('Best Fit Power Law')
    axes.set_ylabel('Response Times (ms)')
    axes.set_xlabel('# Experiences')
    axes.set_xlim(1, max(xdata))
    axes.set_ylim(ylim[0], ylim[1])

def eachBlock_PowerFunTest_andGraph_noTrim(df, xRange, yLow, yHigh, Title):
    import pandas as pd, PythonAnalyses, numpy as np, scipy.stats, PythonAnalyses, matplotlib, matplotlib.pyplot as plt, pylab as pl
    index_matrix = np.zeros((np.size(df['Sub#'].unique()),2))
    amp_matrix = np.zeros((np.size(df['Sub#'].unique()),2))
    ind_part_data_present = np.zeros((np.size(df['Sub#'].unique()),xRange))
    ind_part_data_absent = np.zeros((np.size(df['Sub#'].unique()),xRange))
    CIs = np.zeros(xRange)

    X_data_graph = np.arange(1,xRange+1) #numbers 1-24 for x axis of graph
    counter = 0

    for sub in df['Sub#'].unique():
        sub1_pres = df['RT'][df['Sub#']==sub][df['DistCond']==2]
        sub1_abs = df['RT'][df['Sub#']==sub][df['DistCond']==1]

        x_data1 = df['mergeTrial'][df['Sub#']==sub][df['DistCond']==2]
        x_data2 = df['mergeTrial'][df['Sub#']==sub][df['DistCond']==1]

        PresentFit = PythonAnalyses.Fit_PowerFunction(x_data1, sub1_pres) #fit present data to power func
        index_Present = PresentFit[1]; amp_Present = PresentFit[0];

        AbsentFit = PythonAnalyses.Fit_PowerFunction(x_data2, sub1_abs)
        index_Absent = AbsentFit[1]; amp_Absent = AbsentFit[0]
        
        index_matrix[counter,0] = index_Absent
        index_matrix[counter,1] = index_Present
        amp_matrix[counter,0] = amp_Absent
        amp_matrix[counter,1] = amp_Present
        ind_part_data_absent[counter][:] = amp_Absent*(X_data_graph**index_Absent)
        ind_part_data_present[counter][:] = amp_Present*(X_data_graph**index_Present)
        
        counter+=1

    #print '-----------------------------------------------'
    #print 'c1 ='
    #print index_matrix_Block4_Group2
    #print '-----------------------------------------------'
    #print 'c0 ='
    #print amp_matrix_Block4_Group2

    yFunData_Absent = np.mean(ind_part_data_absent,0)
    yFunData_Present = np.mean(ind_part_data_present,0)

    Tvals = np.ones(xRange)*scipy.stats.ttest_rel(ind_part_data_present, ind_part_data_absent)[0]

    CIs = np.absolute(((np.mean(ind_part_data_present-ind_part_data_absent,axis=0))/ \
        (Tvals)*scipy.stats.t.isf(0.025, len(ind_part_data_present-1))))*(2**(0.05)/2) #loftus & masson + baguley 2011

    fig = plt.figure(figsize=(8,6))
    axes = fig.add_subplot(111)
    axes.plot(X_data_graph,  yFunData_Absent, 'blue', label="Singleton Absent")
    axes.plot(X_data_graph, yFunData_Present, 'green', label="Singleton Present")
    axes.fill_between(X_data_graph, yFunData_Absent+CIs, yFunData_Absent-CIs, color='b', alpha=0.3)
    axes.fill_between(X_data_graph, yFunData_Present+CIs, yFunData_Present-CIs, color='g', alpha=0.3)
    axes.legend()
    axes.set_title(Title)
    axes.set_ylabel('Response Times (ms)')
    axes.set_xlabel('# Experiences')
    axes.set_xlim(0, max(X_data_graph)+1)
    axes.set_ylim(yLow, yHigh)

    normTest1 = scipy.stats.shapiro(index_matrix[:,1]-index_matrix[:,0]) #shapiro-wilk test for normality
    ttest1 = scipy.stats.ttest_rel(index_matrix[:,1], index_matrix[:,0])
    normTest2 = scipy.stats.shapiro(amp_matrix[:,1]-amp_matrix[:,0]) #shapiro-wilk test for normality
    ttest2 = scipy.stats.ttest_rel(amp_matrix[:,1], amp_matrix[:,0])
    print '-----------------------------------------------'
    print 'normality test ='
    print 'index/learning rate/c2/beta'
    print normTest1
    print 'amp/ initial performance/c1/alpha'
    print normTest2
    print '-----------------------------------------------'
    print 'ttest test ='
    print 'index/learning rate/c2/beta'
    print ttest1
    print 'amp/initial performance/c1/alpha'
    print ttest2
    return (index_matrix, amp_matrix)

def eachBlock_PowerFunTest_andGraph_noTrim2(df, xRange, yLow, yHigh, Title):
    import pandas as pd, PythonAnalyses, numpy as np, scipy.stats, PythonAnalyses, matplotlib, matplotlib.pyplot as plt, pylab as pl
    index_matrix = np.zeros((np.size(df['Sub#'].unique()),2))
    amp_matrix = np.zeros((np.size(df['Sub#'].unique()),2))
    ind_part_data_present = np.zeros((np.size(df['Sub#'].unique()),xRange))
    ind_part_data_absent = np.zeros((np.size(df['Sub#'].unique()),xRange))
    CIs = np.zeros(xRange)

    X_data_graph = np.arange(1,xRange+1) #numbers 1-24 for x axis of graph
    counter = 0

    for sub in df['Sub#'].unique():
        sub1_pres = df['RT'][df['Sub#']==sub][df['DistCond']==2]
        sub1_abs = df['RT'][df['Sub#']==sub][df['DistCond']==1]

        x_data1 = df['mergeTrial'][df['Sub#']==sub][df['DistCond']==2]
        x_data2 = df['mergeTrial'][df['Sub#']==sub][df['DistCond']==1]

        PresentFit = PythonAnalyses.Fit_PowerFunction(x_data1, sub1_pres) #fit present data to power func
        index_Present = PresentFit[1]; amp_Present = PresentFit[0];

        AbsentFit = PythonAnalyses.Fit_PowerFunction(x_data2, sub1_abs)
        index_Absent = AbsentFit[1]; amp_Absent = AbsentFit[0]
        
        index_matrix[counter,0] = index_Absent
        index_matrix[counter,1] = index_Present
        amp_matrix[counter,0] = amp_Absent
        amp_matrix[counter,1] = amp_Present
        ind_part_data_absent[counter][:] = amp_Absent*(X_data_graph**index_Absent)
        ind_part_data_present[counter][:] = amp_Present*(X_data_graph**index_Present)
        
        counter+=1

    #print '-----------------------------------------------'
    #print 'c1 ='
    #print index_matrix_Block4_Group2
    #print '-----------------------------------------------'
    #print 'c0 ='
    #print amp_matrix_Block4_Group2

    yFunData_Absent = np.mean(ind_part_data_absent,0)
    yFunData_Present = np.mean(ind_part_data_present,0)

    Tvals = np.ones(xRange)*scipy.stats.ttest_rel(ind_part_data_present, ind_part_data_absent)[0]

    CIs = np.absolute(((np.mean(ind_part_data_present-ind_part_data_absent,axis=0))/ \
        (Tvals)*scipy.stats.t.isf(0.025, len(ind_part_data_present-1))))*(2**(0.05)/2) #loftus & masson + baguley 2011

    df['AvgTrialer'] = np.ceil(df['mergeTrial']/4) #1st half of block. 2nd half. 
    meaner = df.pivot_table(values = 'RT', index='AvgTrialer', columns = ['DistCond'], aggfunc=np.mean)

    fig = plt.figure(figsize=(8,6))
    axes = fig.add_subplot(111)
    axes.plot(np.dot(np.arange(np.size(meaner,0))+1,4)-1.5, meaner[2],'go', markersize=5)
    axes.plot(np.dot(np.arange(np.size(meaner,0))+1,4)-1.5, meaner[1],'bs', markersize=5, markeredgecolor='b', markerfacecolor='None') 
    axes.plot(X_data_graph, yFunData_Present, 'go-', markersize=0.1, markeredgecolor='g', markerfacecolor='g',label="Singleton Present")
    axes.plot(X_data_graph,  yFunData_Absent, 'bs-', markersize=0.1,markeredgecolor='b', markerfacecolor='None',label="Singleton Absent")
    axes.fill_between(X_data_graph, yFunData_Absent+CIs, yFunData_Absent-CIs, color='b', alpha=0.3)
    axes.fill_between(X_data_graph, yFunData_Present+CIs, yFunData_Present-CIs, color='g', alpha=0.3)
    axes.legend(markerscale=50, numpoints=1)
    axes.set_title(Title)
    axes.set_ylabel('Response Times (ms)')
    axes.set_xlabel('# Experiences')
    axes.set_xlim(0, max(X_data_graph)+1)
    axes.set_ylim(yLow, yHigh)

    normTest1 = scipy.stats.shapiro(index_matrix[:,1]-index_matrix[:,0]) #shapiro-wilk test for normality
    ttest1 = scipy.stats.ttest_rel(index_matrix[:,1], index_matrix[:,0])
    normTest2 = scipy.stats.shapiro(amp_matrix[:,1]-amp_matrix[:,0]) #shapiro-wilk test for normality
    ttest2 = scipy.stats.ttest_rel(amp_matrix[:,1], amp_matrix[:,0])
    print '-----------------------------------------------'
    print 'normality test ='
    print 'index/learning rate/c2/beta'
    print normTest1
    print 'amp/ initial performance/c1/alpha'
    print normTest2
    print '-----------------------------------------------'
    print 'ttest test ='
    print 'index/learning rate/c2/beta'
    print ttest1
    print 'amp/initial performance/c1/alpha'
    print ttest2
    return (index_matrix, amp_matrix)

def eachBlock_PowerFunTest_andGraph_noTrim3(df, xRange, yLow, yHigh, Title):
    import pandas as pd, PythonAnalyses, numpy as np, scipy.stats, PythonAnalyses, matplotlib, matplotlib.pyplot as plt, pylab as pl
    index_matrix = np.zeros((np.size(df['Sub#'].unique()),2))
    amp_matrix = np.zeros((np.size(df['Sub#'].unique()),2))
    asmp_matrix = np.zeros((np.size(df['Sub#'].unique()),2))
    ind_part_data_present = np.zeros((np.size(df['Sub#'].unique()),xRange))
    ind_part_data_absent = np.zeros((np.size(df['Sub#'].unique()),xRange))
    CIs = np.zeros(xRange)

    X_data_graph = np.arange(1,xRange+1) #numbers 1-24 for x axis of graph
    counter = 0

    for sub in df['Sub#'].unique():
        sub1_pres = df['RT'][df['Sub#']==sub][df['DistCond']==2]
        sub1_abs = df['RT'][df['Sub#']==sub][df['DistCond']==1]

        x_data1 = df['mergeTrial'][df['Sub#']==sub][df['DistCond']==2]
        x_data2 = df['mergeTrial'][df['Sub#']==sub][df['DistCond']==1]

        PresentFit = PythonAnalyses.Fit_PowerFunction(x_data1, sub1_pres) #fit present data to power func
        index_Present = PresentFit[1]; amp_Present = PresentFit[2]; asmp_Present = PresentFit[3]

        AbsentFit = PythonAnalyses.Fit_PowerFunction(x_data2, sub1_abs)
        index_Absent = AbsentFit[1]; amp_Absent = AbsentFit[2]; asmp_Absent = AbsentFit[3]
        
        index_matrix[counter,0] = index_Absent
        index_matrix[counter,1] = index_Present
        amp_matrix[counter,0] = amp_Absent
        amp_matrix[counter,1] = amp_Present
        asmp_matrix[counter,0] = asmp_Absent
        asmp_matrix[counter,1] = asmp_Present
        ind_part_data_absent[counter][:] = amp_Absent*(X_data_graph**index_Absent)+asmp_Absent
        ind_part_data_present[counter][:] = amp_Present*(X_data_graph**index_Present)+asmp_Present
        
        counter+=1

    #print '-----------------------------------------------'
    #print 'c1 ='
    #print index_matrix_Block4_Group2
    #print '-----------------------------------------------'
    #print 'c0 ='
    #print amp_matrix_Block4_Group2

    yFunData_Absent = np.mean(ind_part_data_absent,0)
    yFunData_Present = np.mean(ind_part_data_present,0)

    Tvals = np.ones(xRange)*scipy.stats.ttest_rel(ind_part_data_present, ind_part_data_absent)[0]

    CIs = np.absolute(((np.mean(ind_part_data_present-ind_part_data_absent,axis=0))/ \
        (Tvals)*scipy.stats.t.isf(0.025, len(ind_part_data_present-1))))*(2**(0.05)/2) #loftus & masson + baguley 2011

    df['AvgTrialer'] = np.ceil(df['mergeTrial']/4) #1st half of block. 2nd half. 
    meaner = df.pivot_table(values = 'RT', index='mergeTrial', columns = ['DistCond'], aggfunc=np.mean)

    fig = plt.figure(figsize=(8,6))
    axes = fig.add_subplot(111)
    axes.plot(np.arange(np.size(meaner,0))+1, meaner[2],'go', markersize=5)
    axes.plot(np.arange(np.size(meaner,0))+1, meaner[1],'bs', markersize=5, markeredgecolor='b', markerfacecolor='None') 
    axes.plot(X_data_graph, yFunData_Present, 'go-', markersize=0.1, markeredgecolor='g', markerfacecolor='g',label="Singleton Present")
    axes.plot(X_data_graph,  yFunData_Absent, 'bs-', markersize=0.1,markeredgecolor='b', markerfacecolor='None',label="Singleton Absent")
    axes.fill_between(X_data_graph, yFunData_Absent+CIs, yFunData_Absent-CIs, color='b', alpha=0.3)
    axes.fill_between(X_data_graph, yFunData_Present+CIs, yFunData_Present-CIs, color='g', alpha=0.3)
    axes.legend(markerscale=50, numpoints=1)
    axes.set_title(Title)
    axes.set_ylabel('Response Times (ms)')
    axes.set_xlabel('# Experiences')
    axes.set_xlim(0, max(X_data_graph)+1)
    axes.set_ylim(yLow, yHigh)

    normTest1 = scipy.stats.shapiro(index_matrix[:,1]-index_matrix[:,0]) #shapiro-wilk test for normality
    ttest1 = scipy.stats.ttest_rel(index_matrix[:,1], index_matrix[:,0])
    normTest2 = scipy.stats.shapiro(amp_matrix[:,1]-amp_matrix[:,0]) #shapiro-wilk test for normality
    ttest2 = scipy.stats.ttest_rel(amp_matrix[:,1], amp_matrix[:,0])
    normTest3 = scipy.stats.shapiro(asmp_matrix[:,1]-asmp_matrix[:,0]) #shapiro-wilk test for normality
    ttest3 = scipy.stats.ttest_rel(asmp_matrix[:,1], asmp_matrix[:,0])
    print '-----------------------------------------------'
    print 'normality test ='
    print 'index/learning rate/c2/beta'
    print normTest1
    print 'amp/ initial performance/c1/alpha'
    print normTest2
    print 'asmp'
    print normTest3
    print '-----------------------------------------------'
    print 'ttest test ='
    print 'index/learning rate/c2/beta'
    print ttest1
    print 'amp/initial performance/c1/alpha'
    print ttest2
    print 'asmp'
    print ttest3
    return (index_matrix, amp_matrix, asmp_matrix)

def rmse_powerlaw(regr,x,y):
    from sklearn import linear_model
    import scipy, numpy as np
    p = regr.predict(scipy.log10(x))
    sum_square_err = sum(np.square(10.0**p - y))
    rmse = np.sqrt(sum_square_err/(len(x)-2))
    return rmse
def rmse_intercept(x,y):
    import numpy as np
    sum_square_err = sum(np.square(np.mean(x)-y))
    rmse = np.sqrt(sum_square_err/len(x)-1)
    return rmse

def rmse_powerlaw2(regr,x,y):
    from sklearn import linear_model
    import scipy, numpy as np
    p = regr.predict(scipy.log10(x))
    sum_square_err = sum(np.square(10.0**p - y))
    rmse = np.sqrt(sum_square_err)
    return rmse
def rmse_intercept2(x,y):
    import numpy as np
    sum_square_err = sum(np.square(np.mean(x)-y))
    rmse = np.sqrt(sum_square_err)
    return rmse

def eachBlock_PowerFunTest_andGraph_noTrim5(df, xRange, yLow, yHigh, Title):
    import pandas as pd, PythonAnalyses, numpy as np, scipy.stats, PythonAnalyses, matplotlib, matplotlib.pyplot as plt, pylab as pl
    from sklearn import linear_model
    index_matrix = np.zeros((np.size(df['Sub#'].unique()),2))
    amp_matrix = np.zeros((np.size(df['Sub#'].unique()),2))
    ind_part_data_present = np.zeros((np.size(df['Sub#'].unique()),xRange))
    ind_part_data_absent = np.zeros((np.size(df['Sub#'].unique()),xRange))
    rmse_pres = np.zeros((np.size(df['Sub#'].unique()),3))
    rmse_abs = np.zeros((np.size(df['Sub#'].unique()),3))
    CIs = np.zeros(xRange)

    X_data_graph = np.arange(1,xRange+1) #numbers 1-24 for x axis of graph
    counter = 0

    for sub in df['Sub#'].unique():
        sub1_pres = df['RT'][df['Sub#']==sub][df['DistCond']==2]
        sub1_abs = df['RT'][df['Sub#']==sub][df['DistCond']==1]

        x_data1 = df['mergeTrial'][df['Sub#']==sub][df['DistCond']==2]
        x_data2 = df['mergeTrial'][df['Sub#']==sub][df['DistCond']==1]

        PresentFit = PythonAnalyses.Fit_PowerFunction(x_data1, sub1_pres) #fit present data to power func
        index_Present = PresentFit[1]; amp_Present = PresentFit[0];

        AbsentFit = PythonAnalyses.Fit_PowerFunction(x_data2, sub1_abs)
        index_Absent = AbsentFit[1]; amp_Absent = AbsentFit[0]

        x_data1 = x_data1.reshape((x_data1.shape[0],-1))
        sub1_pres = sub1_pres.reshape((sub1_pres.shape[0],-1))
        x_data2 = x_data2.reshape((x_data2.shape[0],-1))
        sub1_abs = sub1_abs.reshape((sub1_abs.shape[0],-1))
    
        regr_powerPres = linear_model.LinearRegression()
        regr_powerPres.fit(scipy.log10(x_data1),scipy.log10(sub1_pres))
        rmse_pres[counter,0] = rmse_powerlaw(regr_powerPres,x_data1,sub1_pres)
    
        rmse_pres[counter,1] = rmse_intercept(sub1_pres,sub1_pres)
    
        regr_powerAbs = linear_model.LinearRegression()
        regr_powerAbs.fit(scipy.log10(x_data2),scipy.log10(sub1_abs))
        rmse_abs[counter,0] = rmse_powerlaw(regr_powerAbs,x_data2,sub1_abs)

        rmse_abs[counter,1] = rmse_intercept(sub1_abs,sub1_abs)
        
        index_matrix[counter,0] = index_Absent
        index_matrix[counter,1] = index_Present
        amp_matrix[counter,0] = amp_Absent
        amp_matrix[counter,1] = amp_Present
        ind_part_data_absent[counter][:] = amp_Absent*(X_data_graph**index_Absent)
        ind_part_data_present[counter][:] = amp_Present*(X_data_graph**index_Present)
        
        counter+=1

    rmse_pres[:,2] = rmse_pres[:,0]/rmse_pres[:,1]
    rmse_abs[:,2] = rmse_abs[:,0]/rmse_abs[:,1]

    #yFunData_Absent = np.mean(amp_matrix[:,0])*(X_data_graph**np.mean(index_matrix[:,0]))
    #yFunData_Present = np.mean(amp_matrix[:,1])*(X_data_graph**np.mean(index_matrix[:,1]))

    yFunData_Absent = np.mean(ind_part_data_absent,0)
    yFunData_Present = np.mean(ind_part_data_present,0)

    meaner = df.pivot_table(values = 'RT', index='mergeTrial',columns=['DistCond'], aggfunc=np.mean)
    
    Tvals = np.ones(xRange)*scipy.stats.ttest_rel(ind_part_data_present, ind_part_data_absent)[0]

    CIs = np.absolute(((np.mean(ind_part_data_present-ind_part_data_absent,axis=0))/ \
        (Tvals)*scipy.stats.t.isf(0.025, len(ind_part_data_present-1))))*(2**(0.05)/2) #loftus & masson + baguley 2011

    fig = plt.figure(figsize=(8,6))
    axes = fig.add_subplot(111)
    axes.plot(X_data_graph, meaner[2],'go', markersize=5)
    axes.plot(X_data_graph, meaner[1],'bs', markersize=5, markeredgecolor='b', markerfacecolor='None')
    axes.plot(X_data_graph,  yFunData_Absent, 'bs-', markersize=0.1, markeredgecolor='b', markerfacecolor='None',label="Singleton Absent")
    axes.plot(X_data_graph, yFunData_Present, 'go-', markersize=0.1, markeredgecolor='g', markerfacecolor='g', label="Singleton Present")
    axes.fill_between(X_data_graph, yFunData_Absent+CIs, yFunData_Absent-CIs, color='b', alpha=0.3)
    axes.fill_between(X_data_graph, yFunData_Present+CIs, yFunData_Present-CIs, color='g', alpha=0.3)
    axes.legend(markerscale=50, numpoints=1)
    axes.set_title(Title)
    axes.set_ylabel('Response Times (ms)')
    axes.set_xlabel('# Experiences')
    axes.set_xlim(0, max(X_data_graph)+1)
    axes.set_ylim(yLow, yHigh)

    normTest1 = scipy.stats.shapiro(index_matrix[:,1]-index_matrix[:,0]) #shapiro-wilk test for normality
    ttest1 = scipy.stats.ttest_rel(index_matrix[:,1], index_matrix[:,0])
    normTest2 = scipy.stats.shapiro(amp_matrix[:,1]-amp_matrix[:,0]) #shapiro-wilk test for normality
    ttest2 = scipy.stats.ttest_rel(amp_matrix[:,1], amp_matrix[:,0])
    print 'RT = beta*Trial#^-alpha'
    print '-----------------------------------------------'
    print 'normality test ='
    print 'index/learning rate/c2/alpha'
    print normTest1
    print 'amp/ initial performance/c1/beta'
    print normTest2
    print '-----------------------------------------------'
    print 'ttest test ='
    print 'index/learning rate/c2/alpha'
    print ttest1
    print 'amp/initial performance/c1/beta'
    print ttest2
    return (index_matrix, amp_matrix, rmse_pres, rmse_abs)

def eachBlock_ExpFunTest_andGraph_noTrim(df, xRange, yLow, yHigh, Title, bounders=[(0,2000),(0,100),(0,3000)], guess1=[50.0, 0.1, 850], guess2=[50.0, 0.1, 850]):
    import pandas as pd, PythonAnalyses, numpy as np, scipy.stats, PythonAnalyses, matplotlib, matplotlib.pyplot as plt, pylab as pl
    index_matrix = np.zeros((np.size(df['Sub#'].unique()),2))
    amp_matrix = np.zeros((np.size(df['Sub#'].unique()),2))
    asmp_matrix = np.zeros((np.size(df['Sub#'].unique()),2))
    ind_part_data_present = np.zeros((np.size(df['Sub#'].unique()),xRange))
    ind_part_data_absent = np.zeros((np.size(df['Sub#'].unique()),xRange))
    rmse_pres = np.zeros((np.size(df['Sub#'].unique()),3))
    rmse_abs = np.zeros((np.size(df['Sub#'].unique()),3))
    CIs = np.zeros(xRange)

    X_data_graph = np.arange(1,xRange+1) #numbers 1-24 for x axis of graph
    counter = 0

    for sub in df['Sub#'].unique():
        sub1_pres = df['RT'][df['Sub#']==sub][df['DistCond']==2];
        sub1_abs = df['RT'][df['Sub#']==sub][df['DistCond']==1]; 

        x_data1 = df['mergeTrial'][df['Sub#']==sub][df['DistCond']==2]; 
        x_data2 = df['mergeTrial'][df['Sub#']==sub][df['DistCond']==1];

        PresentFit = PythonAnalyses.Fit_ExpDecFunction(x_data1, sub1_pres, bounders, guess1) #fit present data to power func
        index_Present = PresentFit[1]; amp_Present = PresentFit[2]; asmp_Present = PresentFit[3]

        AbsentFit = PythonAnalyses.Fit_ExpDecFunction(x_data2, sub1_abs, bounders, guess2)
        index_Absent = AbsentFit[1]; amp_Absent = AbsentFit[2]; asmp_Absent = AbsentFit[3]
    
        rmse_pres[counter,0] = PresentFit[0]
        rmse_pres[counter,1] = rmse_intercept(sub1_pres,sub1_pres)
    
        rmse_abs[counter,0] = AbsentFit[0]
        rmse_abs[counter,1] = rmse_intercept(sub1_abs,sub1_abs)
        
        index_matrix[counter,0] = index_Absent
        index_matrix[counter,1] = index_Present
        asmp_matrix[counter,0] = asmp_Absent
        asmp_matrix[counter,1] = asmp_Present
        amp_matrix[counter,0] = amp_Absent
        amp_matrix[counter,1] = amp_Present
        ind_part_data_absent[counter][:] = amp_Absent*np.exp(-index_Absent*X_data_graph)+asmp_Absent
        ind_part_data_present[counter][:] = amp_Present*np.exp(-index_Present*X_data_graph)+asmp_Present
        
        counter+=1

    rmse_pres[:,2] = rmse_pres[:,0]/rmse_pres[:,1]
    rmse_abs[:,2] = rmse_abs[:,0]/rmse_abs[:,1]

    yFunData_Absent = np.mean(ind_part_data_absent,0)
    yFunData_Present = np.mean(ind_part_data_present,0)

    meaner = df.pivot_table(values = 'RT', index='mergeTrial',columns=['DistCond'], aggfunc=np.mean)
    
    Tvals = np.ones(xRange)*scipy.stats.ttest_rel(ind_part_data_present, ind_part_data_absent)[0]

    CIs = np.absolute(((np.mean(ind_part_data_present-ind_part_data_absent,axis=0))/ \
        (Tvals)*scipy.stats.t.isf(0.025, len(ind_part_data_present-1))))*(2**(0.05)/2) #loftus & masson + baguley 2011

    fig = plt.figure(figsize=(8,6))
    axes = fig.add_subplot(111)
    axes.plot(X_data_graph, meaner[2],'go', markersize=5)
    axes.plot(X_data_graph, meaner[1],'bs', markersize=5, markeredgecolor='b', markerfacecolor='None')
    axes.plot(X_data_graph,  yFunData_Absent, 'bs-', markersize=0.1, markeredgecolor='b', markerfacecolor='None',label="Singleton Absent")
    axes.plot(X_data_graph, yFunData_Present, 'go-', markersize=0.1, markeredgecolor='g', markerfacecolor='g', label="Singleton Present")
    axes.fill_between(X_data_graph, yFunData_Absent+CIs, yFunData_Absent-CIs, color='b', alpha=0.3)
    axes.fill_between(X_data_graph, yFunData_Present+CIs, yFunData_Present-CIs, color='g', alpha=0.3)
    axes.legend(markerscale=50, numpoints=1)
    axes.set_title(Title)
    axes.set_ylabel('Response Times (ms)')
    axes.set_xlabel('# Experiences')
    axes.set_xlim(0, max(X_data_graph)+1)
    axes.set_ylim(yLow, yHigh)

    normTest1 = scipy.stats.shapiro(index_matrix[:,1]-index_matrix[:,0]) #shapiro-wilk test for normality
    ttest1 = scipy.stats.ttest_rel(index_matrix[:,1], index_matrix[:,0])
    normTest2 = scipy.stats.shapiro(amp_matrix[:,1]-amp_matrix[:,0]) #shapiro-wilk test for normality
    ttest2 = scipy.stats.ttest_rel(amp_matrix[:,1], amp_matrix[:,0])
    normTest3 = scipy.stats.shapiro(asmp_matrix[:,1]-asmp_matrix[:,0])
    ttest3 = scipy.stats.ttest_rel(asmp_matrix[:,1], asmp_matrix[:,0])
    print 'RT = beta*exp(-alpha*Trial#)+gamma'
    print '-----------------------------------------------'
    print 'normality test ='
    print 'index/learning rate/c2/alpha'
    print normTest1
    print 'amp/ initial performance/c1/beta'
    print normTest2
    print 'asmp/ lower asymptote/c3/gamma'
    print normTest3
    print '-----------------------------------------------'
    print 'ttest test ='
    print 'index/learning rate/c2/alpha'
    print ttest1
    print 'amp/initial performance/c1/beta'
    print ttest2
    print 'asmp/ lower asymptote/c3/gamma'
    print ttest3
    return (index_matrix, amp_matrix, asmp_matrix, rmse_pres, rmse_abs)

def eachBlock_ExpFunTest_noTrim(df, xRange, bounders=[(0,2000),(0,100),(0,3000)], guess1=[50.0, 0.1, 850], guess2=[50.0, 0.1, 850]):
    import pandas as pd, PythonAnalyses, numpy as np, scipy.stats, PythonAnalyses, matplotlib, matplotlib.pyplot as plt, pylab as pl
    index_matrix = np.zeros((np.size(df['Sub#'].unique()),2))
    amp_matrix = np.zeros((np.size(df['Sub#'].unique()),2))
    asmp_matrix = np.zeros((np.size(df['Sub#'].unique()),2))
    rmse_pres = np.zeros((np.size(df['Sub#'].unique()),3))
    rmse_abs = np.zeros((np.size(df['Sub#'].unique()),3))

    X_data_graph = np.arange(1,xRange+1) #numbers 1-24 for x axis of graph
    counter = 0

    for sub in df['Sub#'].unique():
        sub1_pres = df['RT'][df['Sub#']==sub][df['DistCond']==2];
        sub1_abs = df['RT'][df['Sub#']==sub][df['DistCond']==1]; 

        x_data1 = df['mergeTrial'][df['Sub#']==sub][df['DistCond']==2]; 
        x_data2 = df['mergeTrial'][df['Sub#']==sub][df['DistCond']==1];

        PresentFit = PythonAnalyses.Fit_ExpDecFunction(x_data1, sub1_pres, bounders, guess1) #fit present data to power func
        index_Present = PresentFit[1]; amp_Present = PresentFit[2]; asmp_Present = PresentFit[3]

        AbsentFit = PythonAnalyses.Fit_ExpDecFunction(x_data2, sub1_abs, bounders, guess2)
        index_Absent = AbsentFit[1]; amp_Absent = AbsentFit[2]; asmp_Absent = AbsentFit[3]
    
        rmse_pres[counter,0] = PresentFit[0]
        rmse_pres[counter,1] = rmse_intercept(sub1_pres,sub1_pres)
    
        rmse_abs[counter,0] = AbsentFit[0]
        rmse_abs[counter,1] = rmse_intercept(sub1_abs,sub1_abs)
        
        index_matrix[counter,0] = index_Absent
        index_matrix[counter,1] = index_Present
        asmp_matrix[counter,0] = asmp_Present
        asmp_matrix[counter,1] = asmp_Absent
        amp_matrix[counter,0] = amp_Absent
        amp_matrix[counter,1] = amp_Present

        counter+=1

    rmse_pres[:,2] = rmse_pres[:,0]/rmse_pres[:,1]
    rmse_abs[:,2] = rmse_abs[:,0]/rmse_abs[:,1]
