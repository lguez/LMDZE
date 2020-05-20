module paramet_m

  implicit none

  integer, PROTECTED:: iip1, iip2, iip3
  integer, PROTECTED:: jjp1
  integer, PROTECTED:: llmp1, llmp2, llmm1
  integer, PROTECTED:: ip1jm, ip1jmp1, ip1jmi1, ijp1llm, ijmllm

contains

  subroutine paramet

    use dimensions, only: iim, jjm, llm

    !------------------------------------------------------------

    iip1 = iim + 1
    iip2 = iim + 2
    iip3 = iim + 3
    jjp1 = jjm + 1
    llmp1 = llm + 1
    llmp2 = llm + 2
    llmm1 = llm - 1
    ip1jm = (iim + 1) * jjm
    ip1jmp1 = (iim + 1) * (jjm + 1)
    ip1jmi1 = ip1jm - (iim + 1)
    ijp1llm = (iim + 1) * (jjm + 1) * llm
    ijmllm= ip1jm * llm

  end subroutine paramet

end module paramet_m
