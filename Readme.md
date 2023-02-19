# 함수 설명
- 유전자 커널을 생성하기 위한 함수
  - ibs.kernel.intigrate.nobrain(X)
    - IBS커널 계산을 적용해주는 함수
    - 같은 유전자 종류에 포함되는 SNP 데이터 프레임 X를 대입하면 IBS커널로 계산해준다.

  - kernel.matrix.inti.nobrain(X, generange, method = "gaussian", rho = 1)
    - X : 데이터프레임
    - generange : 유전자 데이터가 존재하는 열 넘버
    - method : 커널의 종류 ( 현재 적용되어 있는 함수는 gaussian, ibs )
    - rho : 가우시안 함수를 선택할 경우 조율모수 $\rho$ 적용
    - 결과적으로 리스트 형태의 커널 행렬이 나온다.

- 교호작용 커널을 생성하기 위한 함수
  - kernel.matrix.inti.nobrain(X, generange, brainrange, method = "gaussian", rho = 1)
    - X : 데이터프레임
    - generange : 유전자 데이터가 존재하는 열 넘버
    - brainrange : 관심지역 데이터가 존재하는 열 넘버
    - method : 커널의 종류 ( 현재 적용되어 있는 함수는 gaussian, ibs )
    - rho : 가우시안 함수를 선택할 경우 조율모수 $\rho$ 적용
    - 결과적으로 리스트 형태의 커널 행렬이 나온다.

- 리스트를 하나의 데이터프레임으로 합치는 함수
  - listfordataframe.nobrain(listk)
  - listfordataframe(listk)
    - names가 다르기 때문에 함수를 다르게 정의


