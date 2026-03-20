#!/bin/bash

# ================== 用户可修改部分 ==================
TEMPLATE_FILE="./US_Ne1_soil_acid_P.nml"
OUTPUT_DIR="./namelists/260316CO2/K_REF_CO2"

# ================== 参数检查 ==================
: "${CASE_NAME:?错误: 环境变量 CASE_NAME 未设置}"
: "${FERT:?错误: 环境变量 FERT 未设置}"

# 土壤酸化开关默认值（如果未设置）
: "${USE_SOILACID:=.false.}"
: "${SOILACID_VC:=.false.}"
: "${SOILACID_ER:=.false.}"
: "${SOILACID_NITR:=.false.}"

# 土壤酸化调参默认值（如果未设置）
: "${SOILACID_PH_OPT:=6.5e0}"
: "${SOILACID_PH_SENS:=1.0e0}"
: "${SOILACID_NEU_CONST:=1.55e-9}"
: "${SOILACID_ACID_CONST:=5.0e-4}"
: "${SOILACID_ALK_CONST:=3.3e-7}"
: "${SOILACID_K_REF_CO2:=1.0e-3}"
: "${SOILACID_Z_REF_CO2:=0.5e0}"
: "${SOILACID_EXCHANGE_CALIB:=0.0}"

# 检查模板文件是否存在
if [ ! -f "$TEMPLATE_FILE" ]; then
    echo "错误: 模板文件不存在: $TEMPLATE_FILE"
    exit 1
fi

# 创建输出目录
mkdir -p "$OUTPUT_DIR"

# 生成新的 namelist 文件名
OUTPUT_FILE="${OUTPUT_DIR}/${CASE_NAME}.nml"

# ================== 复制并更新参数 ==================
cp "$TEMPLATE_FILE" "$OUTPUT_FILE"

sed -i -E \
  "s|(DEF_CASE_NAME[[:space:]]*=[[:space:]]*).*$|\1'${CASE_NAME}'|g" \
  "$OUTPUT_FILE"

if [ $? -ne 0 ]; then
    echo "错误: 更新 DEF_CASE_NAME 失败"
    exit 1
fi

sed -i -E \
  "s|(DEF_FERT[[:space:]]*=[[:space:]]*).*$|\1${FERT}|g" \
  "$OUTPUT_FILE"

if [ $? -ne 0 ]; then
    echo "错误: 更新 DEF_FERT 失败"
    exit 1
fi

# ================== 更新土壤酸化开关 ==================
sed -i -E \
  "s|(DEF_USE_SOILACIDIFICATION[[:space:]]*=[[:space:]]*).*$|\1${USE_SOILACID}|g" \
  "$OUTPUT_FILE"

if [ $? -ne 0 ]; then
    echo "错误: 更新 DEF_USE_SOILACIDIFICATION 失败"
    exit 1
fi

sed -i -E \
  "s|(DEF_SOILACID_FACTOR_VC[[:space:]]*=[[:space:]]*).*$|\1${SOILACID_VC}|g" \
  "$OUTPUT_FILE"

if [ $? -ne 0 ]; then
    echo "错误: 更新 DEF_SOILACID_FACTOR_VC 失败"
    exit 1
fi

sed -i -E \
  "s|(DEF_SOILACID_FACTOR_ER[[:space:]]*=[[:space:]]*).*$|\1${SOILACID_ER}|g" \
  "$OUTPUT_FILE"

if [ $? -ne 0 ]; then
    echo "错误: 更新 DEF_SOILACID_FACTOR_ER 失败"
    exit 1
fi

sed -i -E \
  "s|(DEF_SOILACID_FACTOR_NITR[[:space:]]*=[[:space:]]*).*$|\1${SOILACID_NITR}|g" \
  "$OUTPUT_FILE"

if [ $? -ne 0 ]; then
    echo "错误: 更新 DEF_SOILACID_FACTOR_NITR 失败"
    exit 1
fi

# ================== 更新土壤酸化调参 ==================
sed -i -E \
  "s|(DEF_SOILACID_pH_opt[[:space:]]*=[[:space:]]*).*$|\1${SOILACID_PH_OPT}|g" \
  "$OUTPUT_FILE"
if [ $? -ne 0 ]; then echo "错误: 更新 DEF_SOILACID_pH_opt 失败"; exit 1; fi

sed -i -E \
  "s|(DEF_SOILACID_pH_sens[[:space:]]*=[[:space:]]*).*$|\1${SOILACID_PH_SENS}|g" \
  "$OUTPUT_FILE"
if [ $? -ne 0 ]; then echo "错误: 更新 DEF_SOILACID_pH_sens 失败"; exit 1; fi

sed -i -E \
  "s|(DEF_SOILACID_neu_const[[:space:]]*=[[:space:]]*).*$|\1${SOILACID_NEU_CONST}|g" \
  "$OUTPUT_FILE"
if [ $? -ne 0 ]; then echo "错误: 更新 DEF_SOILACID_neu_const 失败"; exit 1; fi

sed -i -E \
  "s|(DEF_SOILACID_acid_const[[:space:]]*=[[:space:]]*).*$|\1${SOILACID_ACID_CONST}|g" \
  "$OUTPUT_FILE"
if [ $? -ne 0 ]; then echo "错误: 更新 DEF_SOILACID_acid_const 失败"; exit 1; fi

sed -i -E \
  "s|(DEF_SOILACID_alk_const[[:space:]]*=[[:space:]]*).*$|\1${SOILACID_ALK_CONST}|g" \
  "$OUTPUT_FILE"
if [ $? -ne 0 ]; then echo "错误: 更新 DEF_SOILACID_alk_const 失败"; exit 1; fi

sed -i -E \
  "s|(DEF_SOILACID_k_ref_co2[[:space:]]*=[[:space:]]*).*$|\1${SOILACID_K_REF_CO2}|g" \
  "$OUTPUT_FILE"
if [ $? -ne 0 ]; then echo "错误: 更新 DEF_SOILACID_k_ref_co2 失败"; exit 1; fi

sed -i -E \
  "s|(DEF_SOILACID_z_ref_co2[[:space:]]*=[[:space:]]*).*$|\1${SOILACID_Z_REF_CO2}|g" \
  "$OUTPUT_FILE"
if [ $? -ne 0 ]; then echo "错误: 更新 DEF_SOILACID_z_ref_co2 失败"; exit 1; fi

sed -i -E \
  "s|(DEF_SOILACID_exchange_calib[[:space:]]*=[[:space:]]*).*$|\1${SOILACID_EXCHANGE_CALIB}|g" \
  "$OUTPUT_FILE"
if [ $? -ne 0 ]; then echo "错误: 更新 DEF_SOILACID_exchange_calib 失败"; exit 1; fi

# ================== 输出确认信息 ==================
echo "========================================"
echo "update_namelist.sh 执行成功"
echo "========================================"
echo "  模板文件       : ${TEMPLATE_FILE}"
echo "  输出文件       : ${OUTPUT_FILE}"
echo "  DEF_CASE_NAME  : ${CASE_NAME}"
echo "  DEF_FERT       : ${FERT}"
echo "----------------------------------------"
echo "  土壤酸化开关:"
echo "    USE_SOILACID      : ${USE_SOILACID}"
echo "    SOILACID_VC       : ${SOILACID_VC}"
echo "    SOILACID_ER       : ${SOILACID_ER}"
echo "    SOILACID_NITR     : ${SOILACID_NITR}"
echo "----------------------------------------"
echo "  土壤酸化调参:"
echo "    PH_OPT            : ${SOILACID_PH_OPT}"
echo "    PH_SENS           : ${SOILACID_PH_SENS}"
echo "    NEU_CONST         : ${SOILACID_NEU_CONST}"
echo "    ACID_CONST        : ${SOILACID_ACID_CONST}"
echo "    ALK_CONST         : ${SOILACID_ALK_CONST}"
echo "    K_REF_CO2         : ${SOILACID_K_REF_CO2}"
echo "    Z_REF_CO2         : ${SOILACID_Z_REF_CO2}"
echo "    EXCHANGE_CALIB    : ${SOILACID_EXCHANGE_CALIB}"
echo "========================================"

# 显示更新后的相关行
echo "更新后的配置内容："
grep -E "DEF_CASE_NAME|DEF_FERT|DEF_USE_SOILACIDIFICATION|DEF_SOILACID" "$OUTPUT_FILE" | sed 's/^/  /'
echo "========================================"

# 导出输出文件路径供后续脚本使用
export NAMELIST_FILE="$OUTPUT_FILE"
echo "已设置环境变量: NAMELIST_FILE=${NAMELIST_FILE}"
