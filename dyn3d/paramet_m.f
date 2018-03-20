module paramet_m

  use dimensions, only: iim, jjm, llm

  implicit none

  private iim, jjm, llm

  integer, PARAMETER:: iip1 = iim + 1, iip2 = iim + 2, iip3 = iim + 3
  integer, PARAMETER:: jjp1 = jjm + 1
  integer, PARAMETER:: llmp1 = llm + 1, llmp2 = llm + 2, llmm1 = llm - 1
  integer, PARAMETER:: ip1jm = (iim + 1) * jjm
  integer, PARAMETER:: ip1jmp1 = (iim + 1) * (jjm + 1)
  integer, PARAMETER:: ip1jmi1 = ip1jm - (iim + 1)
  integer, PARAMETER:: ijp1llm = (iim + 1) * (jjm + 1) * llm
  integer, PARAMETER:: ijmllm= ip1jm * llm 

end module paramet_m
