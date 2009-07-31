test_ozonecm.o : phyetat0_test.o ozonecm.o comvert.o dimens_m_test.o 
comvert.o : dimens_m_test.o 
ozonecm.o : phyetat0_test.o dimphy_test.o dimens_m_test.o 
phyetat0_test.o : dimphy_test.o 
dimphy_test.o : dimens_m_test.o 
