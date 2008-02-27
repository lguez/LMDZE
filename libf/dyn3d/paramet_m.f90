module paramet_m

  use dimens_m, only: iim, jjm, llm, ndm

  implicit none

  private iim, jjm, llm, ndm

  integer, PARAMETER:: iip1 = iim + 1, iip2 = iim + 2, iip3 = iim + 3
  integer, PARAMETER:: jjp1 = jjm + 1
  integer, PARAMETER:: llmp1 = llm + 1,  llmp2 = llm + 2, llmm1 = llm - 1

  integer, PARAMETER:: kftd  = iim/2 -ndm 
  integer, PARAMETER:: ip1jm  = (iim + 1) * jjm
  integer, PARAMETER:: ip1jmp1= (iim + 1) * (jjm + 1) 
  integer, PARAMETER:: ip1jmi1= ip1jm - (iim + 1) 
  integer, PARAMETER:: ijp1llm= (iim + 1) * (jjm + 1) * llm
  integer, PARAMETER:: ijmllm= ip1jm * llm 
  integer, PARAMETER:: mvar= (iim + 1) * (jjm + 1)*( 2*llm+1) + ijmllm 
  integer, PARAMETER:: jcfil=jjm/2+5, jcfllm=jcfil*llm 

end module paramet_m
