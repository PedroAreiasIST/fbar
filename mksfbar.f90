
SUBROUTINE mksfbar(v,c,cb,cs,dcsdc,dcsdc2)
USE SMSUtility
IMPLICIT NONE
DOUBLE PRECISION v(1716),c(6),cb(6),cs(6),dcsdc(6,12),dcsdc2(6,12,12)
v(1710)=(-2d0)*c(6)
v(1706)=(-2d0)*cb(6)
v(1699)=cb(5)**2
v(1698)=cb(4)**2
v(1697)=cb(2)*cb(3)-cb(6)**2
v(1696)=cb(1)*v(1706)
v(1695)=2d0*cb(5)
v(1694)=2d0*cb(4)
v(1700)=cb(5)*v(1694)
v(1693)=cb(6)*v(1694)-cb(2)*v(1695)
v(1692)=cb(1)*cb(3)-v(1699)
v(1691)=cb(1)*cb(2)-v(1698)
v(1689)=c(5)**2
v(1688)=c(4)**2
v(1687)=c(2)*c(3)-c(6)**2
v(1686)=c(1)*v(1710)
v(1685)=2d0*c(5)
v(1684)=2d0*c(4)
v(1690)=c(5)*v(1684)
v(1683)=c(6)*v(1684)-c(2)*v(1685)
v(1682)=c(1)*c(3)-v(1689)
v(1681)=c(1)*c(2)-v(1688)
v(1058)=-(c(3)*v(1684))+c(6)*v(1685)
v(1060)=v(1686)+v(1690)
v(964)=c(1)*v(1687)-c(3)*v(1688)-c(2)*v(1689)+c(6)*v(1690)
v(1061)=1d0/v(964)**2
v(1067)=-(v(1060)*v(1061))
v(1066)=-(v(1061)*v(1683))
v(1065)=-(v(1058)*v(1061))
v(1064)=-(v(1061)*v(1681))
v(1063)=-(v(1061)*v(1682))
v(1062)=-(v(1061)*v(1687))
v(1071)=-(cb(3)*v(1694))+cb(6)*v(1695)
v(1073)=v(1696)+v(1700)
v(965)=cb(1)*v(1697)-cb(3)*v(1698)-cb(2)*v(1699)+cb(6)*v(1700)
v(1703)=(2d0/3d0)*v(965)
v(1171)=(-1d0/3d0)*(v(1061)*v(965))
v(962)=v(965)/v(964)
v(1087)=1d0/v(962)**0.16666666666666669d1
v(1702)=-(v(1087)*v(1703))
v(1701)=(-2d0/9d0)*v(1061)*v(1087)
v(1109)=v(1693)*v(1701)
v(1108)=v(1071)*v(1701)
v(1107)=v(1691)*v(1701)
v(1155)=v(1107)*v(1691)
v(1132)=v(1107)*v(1693)
v(1119)=v(1073)*v(1107)
v(1106)=v(1692)*v(1701)
v(1163)=v(1106)*v(1692)
v(1143)=v(1071)*v(1106)
v(1118)=v(1073)*v(1106)
v(1105)=v(1697)*v(1701)
v(1170)=v(1105)*v(1697)
v(1142)=v(1071)*v(1105)
v(1129)=v(1105)*v(1693)
v(1092)=v(1067)*v(1702)
v(1091)=v(1066)*v(1702)
v(1090)=v(1065)*v(1702)
v(1089)=v(1064)*v(1702)
v(1088)=v(1063)*v(1702)
v(1086)=v(1062)*v(1702)
v(973)=1d0/v(962)**0.6666666666666666d0
v(1173)=(v(1703)*v(973))/v(964)**3
v(1177)=v(1091)*v(1171)+v(1173)*v(1683)
v(1176)=v(1090)*v(1171)+v(1058)*v(1173)
v(1175)=v(1089)*v(1171)+v(1173)*v(1681)
v(1206)=v(1175)*v(1681)
v(1195)=v(1175)*v(1683)
v(1188)=v(1060)*v(1175)
v(1174)=v(1088)*v(1171)+v(1173)*v(1682)
v(1208)=v(1174)*v(1682)
v(1200)=v(1058)*v(1174)
v(1187)=v(1060)*v(1174)
v(1172)=v(1086)*v(1171)+v(1173)*v(1687)
v(1209)=v(1172)*v(1687)
v(1199)=v(1058)*v(1172)
v(1192)=v(1172)*v(1683)
v(1104)=(v(1092)/v(964)+v(1067)*v(973))/3d0
v(1169)=v(1104)*v(1697)
v(1161)=v(1104)*v(1692)
v(1152)=v(1104)*v(1691)
v(1141)=v(1071)*v(1104)
v(1128)=v(1104)*v(1693)
v(1116)=v(1073)*v(1104)
v(1103)=(v(1091)/v(964)+v(1066)*v(973))/3d0
v(1168)=v(1103)*v(1697)
v(1160)=v(1103)*v(1692)
v(1151)=v(1103)*v(1691)
v(1140)=v(1071)*v(1103)
v(1127)=v(1103)*v(1693)
v(1115)=v(1073)*v(1103)
v(1102)=(v(1090)/v(964)+v(1065)*v(973))/3d0
v(1167)=v(1102)*v(1697)
v(1159)=v(1102)*v(1692)
v(1150)=v(1102)*v(1691)
v(1139)=v(1071)*v(1102)
v(1126)=v(1102)*v(1693)
v(1114)=v(1073)*v(1102)
v(1101)=(v(1089)/v(964)+v(1064)*v(973))/3d0
v(1166)=v(1101)*v(1697)
v(1158)=v(1101)*v(1692)
v(1149)=v(1101)*v(1691)
v(1138)=v(1071)*v(1101)
v(1125)=v(1101)*v(1693)
v(1113)=v(1073)*v(1101)
v(1100)=(v(1088)/v(964)+v(1063)*v(973))/3d0
v(1165)=v(1100)*v(1697)
v(1157)=v(1100)*v(1692)
v(1148)=v(1100)*v(1691)
v(1137)=v(1071)*v(1100)
v(1124)=v(1100)*v(1693)
v(1112)=v(1073)*v(1100)
v(1099)=(v(1086)/v(964)+v(1062)*v(973))/3d0
v(1164)=v(1099)*v(1697)
v(1156)=v(1099)*v(1692)
v(1147)=v(1099)*v(1691)
v(1136)=v(1071)*v(1099)
v(1123)=v(1099)*v(1693)
v(1111)=v(1073)*v(1099)
v(977)=v(973)/(3d0*v(964))
v(1707)=cb(1)*v(977)
v(1705)=cb(2)*v(977)
v(1704)=cb(3)*v(977)
v(1162)=v(1105)*v(1692)+v(1704)
v(1154)=v(1106)*v(1691)+v(1707)
v(1153)=v(1105)*v(1691)+v(1705)
v(1146)=v(1071)*v(1108)-2d0*v(1704)
v(1144)=v(1694)*v(977)
v(1145)=v(1071)*v(1107)-v(1144)
v(1135)=v(1109)*v(1693)-2d0*v(1705)
v(1133)=v(1706)*v(977)
v(1134)=-v(1133)+v(1108)*v(1693)
v(1130)=v(1695)*v(977)
v(1131)=-v(1130)+v(1106)*v(1693)
v(1122)=(v(1073)*v(1073))*v(1701)-2d0*v(1707)
v(1121)=v(1073)*v(1109)+v(1144)
v(1120)=v(1073)*v(1108)+v(1130)
v(1117)=v(1073)*v(1105)+v(1133)
v(982)=v(1073)*v(977)
v(981)=v(1693)*v(977)
v(980)=v(1071)*v(977)
v(979)=v(1691)*v(977)
v(978)=v(1692)*v(977)
v(976)=v(1697)*v(977)
v(967)=v(1171)*v(973)
v(1711)=c(1)*v(967)
v(1709)=c(2)*v(967)
v(1708)=c(3)*v(967)
v(1207)=v(1172)*v(1682)+v(1708)
v(1205)=v(1174)*v(1681)+v(1711)
v(1204)=v(1172)*v(1681)+v(1709)
v(1203)=v(1058)*v(1176)-2d0*v(1708)
v(1201)=v(1684)*v(967)
v(1202)=v(1058)*v(1175)-v(1201)
v(1198)=v(1177)*v(1683)-2d0*v(1709)
v(1196)=v(1710)*v(967)
v(1197)=-v(1196)+v(1176)*v(1683)
v(1193)=v(1685)*v(967)
v(1194)=-v(1193)+v(1174)*v(1683)
v(1191)=v(1060)*(v(1092)*v(1171)+v(1060)*v(1173))-2d0*v(1711)
v(1190)=v(1060)*v(1177)+v(1201)
v(1189)=v(1060)*v(1176)+v(1193)
v(1186)=v(1060)*v(1172)+v(1196)
v(972)=v(1060)*v(967)
v(971)=v(1683)*v(967)
v(970)=v(1058)*v(967)
v(969)=v(1681)*v(967)
v(968)=v(1682)*v(967)
v(966)=v(1687)*v(967)
v(956)=v(962)**0.3333333333333333d0
v(1211)=c(1)*v(1207)+v(968)
v(1212)=c(1)*v(1204)+v(969)
v(1213)=c(1)*v(1199)+v(970)
v(1214)=c(1)*v(1192)+v(971)
v(1215)=c(1)*v(1186)+v(972)
v(1216)=c(1)*v(1164)+v(976)
v(1217)=c(1)*v(1156)+v(978)
v(1218)=c(1)*v(1147)+v(979)
v(1219)=c(1)*v(1136)+v(980)
v(1220)=c(1)*v(1123)+v(981)
v(1221)=c(1)*v(1111)+v(982)
v(1223)=c(2)*v(1207)+v(966)
v(1224)=c(2)*v(1204)
v(1225)=c(2)*v(1199)
v(1226)=c(2)*v(1192)
v(1227)=c(2)*v(1186)
v(1228)=c(2)*v(1164)
v(1229)=c(2)*v(1156)
v(1230)=c(2)*v(1147)
v(1231)=c(2)*v(1136)
v(1232)=c(2)*v(1123)
v(1233)=c(2)*v(1111)
v(1235)=c(3)*v(1207)
v(1236)=c(3)*v(1204)+v(966)
v(1237)=c(3)*v(1199)
v(1238)=c(3)*v(1192)
v(1239)=c(3)*v(1186)
v(1240)=c(3)*v(1164)
v(1241)=c(3)*v(1156)
v(1242)=c(3)*v(1147)
v(1243)=c(3)*v(1136)
v(1244)=c(3)*v(1123)
v(1245)=c(3)*v(1111)
v(1247)=c(4)*v(1207)
v(1248)=c(4)*v(1204)
v(1249)=c(4)*v(1199)+v(966)
v(1250)=c(4)*v(1192)
v(1251)=c(4)*v(1186)
v(1252)=c(4)*v(1164)
v(1253)=c(4)*v(1156)
v(1254)=c(4)*v(1147)
v(1255)=c(4)*v(1136)
v(1256)=c(4)*v(1123)
v(1257)=c(4)*v(1111)
v(1259)=c(5)*v(1207)
v(1260)=c(5)*v(1204)
v(1261)=c(5)*v(1199)
v(1262)=c(5)*v(1192)+v(966)
v(1263)=c(5)*v(1186)
v(1264)=c(5)*v(1164)
v(1265)=c(5)*v(1156)
v(1266)=c(5)*v(1147)
v(1267)=c(5)*v(1136)
v(1268)=c(5)*v(1123)
v(1269)=c(5)*v(1111)
v(1271)=c(6)*v(1207)
v(1272)=c(6)*v(1204)
v(1273)=c(6)*v(1199)
v(1274)=c(6)*v(1192)
v(1275)=c(6)*v(1186)+v(966)
v(1276)=c(6)*v(1164)
v(1277)=c(6)*v(1156)
v(1278)=c(6)*v(1147)
v(1279)=c(6)*v(1136)
v(1280)=c(6)*v(1123)
v(1281)=c(6)*v(1111)
v(1283)=c(1)*v(1205)
v(1284)=c(1)*v(1200)
v(1285)=c(1)*v(1194)
v(1286)=c(1)*v(1187)
v(1287)=c(1)*v(1165)
v(1288)=c(1)*v(1157)
v(1289)=c(1)*v(1148)
v(1290)=c(1)*v(1137)
v(1291)=c(1)*v(1124)
v(1292)=c(1)*v(1112)
v(1294)=c(2)*v(1205)+v(969)
v(1295)=c(2)*v(1200)+v(970)
v(1296)=c(2)*v(1194)+v(971)
v(1297)=c(2)*v(1187)+v(972)
v(1298)=c(2)*v(1165)+v(976)
v(1299)=c(2)*v(1157)+v(978)
v(1300)=c(2)*v(1148)+v(979)
v(1301)=c(2)*v(1137)+v(980)
v(1302)=c(2)*v(1124)+v(981)
v(1303)=c(2)*v(1112)+v(982)
v(1305)=c(3)*v(1205)+v(968)
v(1306)=c(3)*v(1200)
v(1307)=c(3)*v(1194)
v(1308)=c(3)*v(1187)
v(1309)=c(3)*v(1165)
v(1310)=c(3)*v(1157)
v(1311)=c(3)*v(1148)
v(1312)=c(3)*v(1137)
v(1313)=c(3)*v(1124)
v(1314)=c(3)*v(1112)
v(1316)=c(4)*v(1205)
v(1317)=c(4)*v(1200)+v(968)
v(1318)=c(4)*v(1194)
v(1319)=c(4)*v(1187)
v(1320)=c(4)*v(1165)
v(1321)=c(4)*v(1157)
v(1322)=c(4)*v(1148)
v(1323)=c(4)*v(1137)
v(1324)=c(4)*v(1124)
v(1325)=c(4)*v(1112)
v(1327)=c(5)*v(1205)
v(1328)=c(5)*v(1200)
v(1329)=c(5)*v(1194)+v(968)
v(1330)=c(5)*v(1187)
v(1331)=c(5)*v(1165)
v(1332)=c(5)*v(1157)
v(1333)=c(5)*v(1148)
v(1334)=c(5)*v(1137)
v(1335)=c(5)*v(1124)
v(1336)=c(5)*v(1112)
v(1338)=c(6)*v(1205)
v(1339)=c(6)*v(1200)
v(1340)=c(6)*v(1194)
v(1341)=c(6)*v(1187)+v(968)
v(1342)=c(6)*v(1165)
v(1343)=c(6)*v(1157)
v(1344)=c(6)*v(1148)
v(1345)=c(6)*v(1137)
v(1346)=c(6)*v(1124)
v(1347)=c(6)*v(1112)
v(1349)=c(1)*v(1202)
v(1350)=c(1)*v(1195)
v(1351)=c(1)*v(1188)
v(1352)=c(1)*v(1166)
v(1353)=c(1)*v(1158)
v(1354)=c(1)*v(1149)
v(1355)=c(1)*v(1138)
v(1356)=c(1)*v(1125)
v(1357)=c(1)*v(1113)
v(1359)=c(2)*v(1202)
v(1360)=c(2)*v(1195)
v(1361)=c(2)*v(1188)
v(1362)=c(2)*v(1166)
v(1363)=c(2)*v(1158)
v(1364)=c(2)*v(1149)
v(1365)=c(2)*v(1138)
v(1366)=c(2)*v(1125)
v(1367)=c(2)*v(1113)
v(1369)=c(3)*v(1202)+v(970)
v(1370)=c(3)*v(1195)+v(971)
v(1371)=c(3)*v(1188)+v(972)
v(1372)=c(3)*v(1166)+v(976)
v(1373)=c(3)*v(1158)+v(978)
v(1374)=c(3)*v(1149)+v(979)
v(1375)=c(3)*v(1138)+v(980)
v(1376)=c(3)*v(1125)+v(981)
v(1377)=c(3)*v(1113)+v(982)
v(1379)=c(4)*v(1202)+v(969)
v(1380)=c(4)*v(1195)
v(1381)=c(4)*v(1188)
v(1382)=c(4)*v(1166)
v(1383)=c(4)*v(1158)
v(1384)=c(4)*v(1149)
v(1385)=c(4)*v(1138)
v(1386)=c(4)*v(1125)
v(1387)=c(4)*v(1113)
v(1389)=c(5)*v(1202)
v(1390)=c(5)*v(1195)+v(969)
v(1391)=c(5)*v(1188)
v(1392)=c(5)*v(1166)
v(1393)=c(5)*v(1158)
v(1394)=c(5)*v(1149)
v(1395)=c(5)*v(1138)
v(1396)=c(5)*v(1125)
v(1397)=c(5)*v(1113)
v(1399)=c(6)*v(1202)
v(1400)=c(6)*v(1195)
v(1401)=c(6)*v(1188)+v(969)
v(1402)=c(6)*v(1166)
v(1403)=c(6)*v(1158)
v(1404)=c(6)*v(1149)
v(1405)=c(6)*v(1138)
v(1406)=c(6)*v(1125)
v(1407)=c(6)*v(1113)
v(1409)=c(1)*v(1197)
v(1410)=c(1)*v(1189)
v(1411)=c(1)*v(1167)
v(1412)=c(1)*v(1159)
v(1413)=c(1)*v(1150)
v(1414)=c(1)*v(1139)
v(1415)=c(1)*v(1126)
v(1416)=c(1)*v(1114)
v(1418)=c(2)*v(1197)
v(1419)=c(2)*v(1189)
v(1420)=c(2)*v(1167)
v(1421)=c(2)*v(1159)
v(1422)=c(2)*v(1150)
v(1423)=c(2)*v(1139)
v(1424)=c(2)*v(1126)
v(1425)=c(2)*v(1114)
v(1427)=c(3)*v(1197)
v(1428)=c(3)*v(1189)
v(1429)=c(3)*v(1167)
v(1430)=c(3)*v(1159)
v(1431)=c(3)*v(1150)
v(1432)=c(3)*v(1139)
v(1433)=c(3)*v(1126)
v(1434)=c(3)*v(1114)
v(1436)=c(4)*v(1197)+v(971)
v(1437)=c(4)*v(1189)+v(972)
v(1438)=c(4)*v(1167)+v(976)
v(1439)=c(4)*v(1159)+v(978)
v(1440)=c(4)*v(1150)+v(979)
v(1441)=c(4)*v(1139)+v(980)
v(1442)=c(4)*v(1126)+v(981)
v(1443)=c(4)*v(1114)+v(982)
v(1445)=c(5)*v(1197)+v(970)
v(1446)=c(5)*v(1189)
v(1447)=c(5)*v(1167)
v(1448)=c(5)*v(1159)
v(1449)=c(5)*v(1150)
v(1450)=c(5)*v(1139)
v(1451)=c(5)*v(1126)
v(1452)=c(5)*v(1114)
v(1454)=c(6)*v(1197)
v(1455)=c(6)*v(1189)+v(970)
v(1456)=c(6)*v(1167)
v(1457)=c(6)*v(1159)
v(1458)=c(6)*v(1150)
v(1459)=c(6)*v(1139)
v(1460)=c(6)*v(1126)
v(1461)=c(6)*v(1114)
v(1463)=c(1)*v(1190)
v(1464)=c(1)*v(1168)
v(1465)=c(1)*v(1160)
v(1466)=c(1)*v(1151)
v(1467)=c(1)*v(1140)
v(1468)=c(1)*v(1127)
v(1469)=c(1)*v(1115)
v(1471)=c(2)*v(1190)
v(1472)=c(2)*v(1168)
v(1473)=c(2)*v(1160)
v(1474)=c(2)*v(1151)
v(1475)=c(2)*v(1140)
v(1476)=c(2)*v(1127)
v(1477)=c(2)*v(1115)
v(1479)=c(3)*v(1190)
v(1480)=c(3)*v(1168)
v(1481)=c(3)*v(1160)
v(1482)=c(3)*v(1151)
v(1483)=c(3)*v(1140)
v(1484)=c(3)*v(1127)
v(1485)=c(3)*v(1115)
v(1487)=c(4)*v(1190)
v(1488)=c(4)*v(1168)
v(1489)=c(4)*v(1160)
v(1490)=c(4)*v(1151)
v(1491)=c(4)*v(1140)
v(1492)=c(4)*v(1127)
v(1493)=c(4)*v(1115)
v(1495)=c(5)*v(1190)+v(972)
v(1496)=c(5)*v(1168)+v(976)
v(1497)=c(5)*v(1160)+v(978)
v(1498)=c(5)*v(1151)+v(979)
v(1499)=c(5)*v(1140)+v(980)
v(1500)=c(5)*v(1127)+v(981)
v(1501)=c(5)*v(1115)+v(982)
v(1503)=c(6)*v(1190)+v(971)
v(1504)=c(6)*v(1168)
v(1505)=c(6)*v(1160)
v(1506)=c(6)*v(1151)
v(1507)=c(6)*v(1140)
v(1508)=c(6)*v(1127)
v(1509)=c(6)*v(1115)
v(1511)=c(1)*v(1169)
v(1512)=c(1)*v(1161)
v(1513)=c(1)*v(1152)
v(1514)=c(1)*v(1141)
v(1515)=c(1)*v(1128)
v(1516)=c(1)*v(1116)
v(1518)=c(2)*v(1169)
v(1519)=c(2)*v(1161)
v(1520)=c(2)*v(1152)
v(1521)=c(2)*v(1141)
v(1522)=c(2)*v(1128)
v(1523)=c(2)*v(1116)
v(1525)=c(3)*v(1169)
v(1526)=c(3)*v(1161)
v(1527)=c(3)*v(1152)
v(1528)=c(3)*v(1141)
v(1529)=c(3)*v(1128)
v(1530)=c(3)*v(1116)
v(1532)=c(4)*v(1169)
v(1533)=c(4)*v(1161)
v(1534)=c(4)*v(1152)
v(1535)=c(4)*v(1141)
v(1536)=c(4)*v(1128)
v(1537)=c(4)*v(1116)
v(1539)=c(5)*v(1169)
v(1540)=c(5)*v(1161)
v(1541)=c(5)*v(1152)
v(1542)=c(5)*v(1141)
v(1543)=c(5)*v(1128)
v(1544)=c(5)*v(1116)
v(1546)=c(6)*v(1169)+v(976)
v(1547)=c(6)*v(1161)+v(978)
v(1548)=c(6)*v(1152)+v(979)
v(1549)=c(6)*v(1141)+v(980)
v(1550)=c(6)*v(1128)+v(981)
v(1551)=c(6)*v(1116)+v(982)
v(1553)=c(1)*v(1162)
v(1554)=c(1)*v(1153)
v(1555)=c(1)*v(1142)
v(1556)=c(1)*v(1129)
v(1557)=c(1)*v(1117)
v(1559)=c(2)*v(1162)
v(1560)=c(2)*v(1153)
v(1561)=c(2)*v(1142)
v(1562)=c(2)*v(1129)
v(1563)=c(2)*v(1117)
v(1565)=c(3)*v(1162)
v(1566)=c(3)*v(1153)
v(1567)=c(3)*v(1142)
v(1568)=c(3)*v(1129)
v(1569)=c(3)*v(1117)
v(1571)=c(4)*v(1162)
v(1572)=c(4)*v(1153)
v(1573)=c(4)*v(1142)
v(1574)=c(4)*v(1129)
v(1575)=c(4)*v(1117)
v(1577)=c(5)*v(1162)
v(1578)=c(5)*v(1153)
v(1579)=c(5)*v(1142)
v(1580)=c(5)*v(1129)
v(1581)=c(5)*v(1117)
v(1583)=c(6)*v(1162)
v(1584)=c(6)*v(1153)
v(1585)=c(6)*v(1142)
v(1586)=c(6)*v(1129)
v(1587)=c(6)*v(1117)
v(1589)=c(1)*v(1154)
v(1590)=c(1)*v(1143)
v(1591)=c(1)*v(1131)
v(1592)=c(1)*v(1118)
v(1594)=c(2)*v(1154)
v(1595)=c(2)*v(1143)
v(1596)=c(2)*v(1131)
v(1597)=c(2)*v(1118)
v(1599)=c(3)*v(1154)
v(1600)=c(3)*v(1143)
v(1601)=c(3)*v(1131)
v(1602)=c(3)*v(1118)
v(1604)=c(4)*v(1154)
v(1605)=c(4)*v(1143)
v(1606)=c(4)*v(1131)
v(1607)=c(4)*v(1118)
v(1609)=c(5)*v(1154)
v(1610)=c(5)*v(1143)
v(1611)=c(5)*v(1131)
v(1612)=c(5)*v(1118)
v(1614)=c(6)*v(1154)
v(1615)=c(6)*v(1143)
v(1616)=c(6)*v(1131)
v(1617)=c(6)*v(1118)
v(1619)=c(1)*v(1145)
v(1620)=c(1)*v(1132)
v(1621)=c(1)*v(1119)
v(1623)=c(2)*v(1145)
v(1624)=c(2)*v(1132)
v(1625)=c(2)*v(1119)
v(1627)=c(3)*v(1145)
v(1628)=c(3)*v(1132)
v(1629)=c(3)*v(1119)
v(1631)=c(4)*v(1145)
v(1632)=c(4)*v(1132)
v(1633)=c(4)*v(1119)
v(1635)=c(5)*v(1145)
v(1636)=c(5)*v(1132)
v(1637)=c(5)*v(1119)
v(1639)=c(6)*v(1145)
v(1640)=c(6)*v(1132)
v(1641)=c(6)*v(1119)
v(1643)=c(1)*v(1134)
v(1644)=c(1)*v(1120)
v(1646)=c(2)*v(1134)
v(1647)=c(2)*v(1120)
v(1649)=c(3)*v(1134)
v(1650)=c(3)*v(1120)
v(1652)=c(4)*v(1134)
v(1653)=c(4)*v(1120)
v(1655)=c(5)*v(1134)
v(1656)=c(5)*v(1120)
v(1658)=c(6)*v(1134)
v(1659)=c(6)*v(1120)
v(1661)=c(1)*v(1121)
v(1663)=c(2)*v(1121)
v(1665)=c(3)*v(1121)
v(1667)=c(4)*v(1121)
v(1669)=c(5)*v(1121)
v(1671)=c(6)*v(1121)
cs(1)=c(1)*v(956)
cs(2)=c(2)*v(956)
cs(3)=c(3)*v(956)
cs(4)=c(4)*v(956)
cs(5)=c(5)*v(956)
cs(6)=c(6)*v(956)
dcsdc(1,1)=v(956)+c(1)*v(966)
dcsdc(1,2)=c(1)*v(968)
dcsdc(1,3)=c(1)*v(969)
dcsdc(1,4)=c(1)*v(970)
dcsdc(1,5)=c(1)*v(971)
dcsdc(1,6)=c(1)*v(972)
dcsdc(1,7)=c(1)*v(976)
dcsdc(1,8)=c(1)*v(978)
dcsdc(1,9)=c(1)*v(979)
dcsdc(1,10)=c(1)*v(980)
dcsdc(1,11)=c(1)*v(981)
dcsdc(1,12)=c(1)*v(982)
dcsdc(2,1)=c(2)*v(966)
dcsdc(2,2)=v(956)+c(2)*v(968)
dcsdc(2,3)=c(2)*v(969)
dcsdc(2,4)=c(2)*v(970)
dcsdc(2,5)=c(2)*v(971)
dcsdc(2,6)=c(2)*v(972)
dcsdc(2,7)=c(2)*v(976)
dcsdc(2,8)=c(2)*v(978)
dcsdc(2,9)=c(2)*v(979)
dcsdc(2,10)=c(2)*v(980)
dcsdc(2,11)=c(2)*v(981)
dcsdc(2,12)=c(2)*v(982)
dcsdc(3,1)=c(3)*v(966)
dcsdc(3,2)=c(3)*v(968)
dcsdc(3,3)=v(956)+c(3)*v(969)
dcsdc(3,4)=c(3)*v(970)
dcsdc(3,5)=c(3)*v(971)
dcsdc(3,6)=c(3)*v(972)
dcsdc(3,7)=c(3)*v(976)
dcsdc(3,8)=c(3)*v(978)
dcsdc(3,9)=c(3)*v(979)
dcsdc(3,10)=c(3)*v(980)
dcsdc(3,11)=c(3)*v(981)
dcsdc(3,12)=c(3)*v(982)
dcsdc(4,1)=c(4)*v(966)
dcsdc(4,2)=c(4)*v(968)
dcsdc(4,3)=c(4)*v(969)
dcsdc(4,4)=v(956)+c(4)*v(970)
dcsdc(4,5)=c(4)*v(971)
dcsdc(4,6)=c(4)*v(972)
dcsdc(4,7)=c(4)*v(976)
dcsdc(4,8)=c(4)*v(978)
dcsdc(4,9)=c(4)*v(979)
dcsdc(4,10)=c(4)*v(980)
dcsdc(4,11)=c(4)*v(981)
dcsdc(4,12)=c(4)*v(982)
dcsdc(5,1)=c(5)*v(966)
dcsdc(5,2)=c(5)*v(968)
dcsdc(5,3)=c(5)*v(969)
dcsdc(5,4)=c(5)*v(970)
dcsdc(5,5)=v(956)+c(5)*v(971)
dcsdc(5,6)=c(5)*v(972)
dcsdc(5,7)=c(5)*v(976)
dcsdc(5,8)=c(5)*v(978)
dcsdc(5,9)=c(5)*v(979)
dcsdc(5,10)=c(5)*v(980)
dcsdc(5,11)=c(5)*v(981)
dcsdc(5,12)=c(5)*v(982)
dcsdc(6,1)=c(6)*v(966)
dcsdc(6,2)=c(6)*v(968)
dcsdc(6,3)=c(6)*v(969)
dcsdc(6,4)=c(6)*v(970)
dcsdc(6,5)=c(6)*v(971)
dcsdc(6,6)=v(956)+c(6)*v(972)
dcsdc(6,7)=c(6)*v(976)
dcsdc(6,8)=c(6)*v(978)
dcsdc(6,9)=c(6)*v(979)
dcsdc(6,10)=c(6)*v(980)
dcsdc(6,11)=c(6)*v(981)
dcsdc(6,12)=c(6)*v(982)
dcsdc2(1,1,1)=c(1)*v(1209)+2d0*v(966)
dcsdc2(1,1,2)=v(1211)
dcsdc2(1,1,3)=v(1212)
dcsdc2(1,1,4)=v(1213)
dcsdc2(1,1,5)=v(1214)
dcsdc2(1,1,6)=v(1215)
dcsdc2(1,1,7)=v(1216)
dcsdc2(1,1,8)=v(1217)
dcsdc2(1,1,9)=v(1218)
dcsdc2(1,1,10)=v(1219)
dcsdc2(1,1,11)=v(1220)
dcsdc2(1,1,12)=v(1221)
dcsdc2(1,2,1)=v(1211)
dcsdc2(1,2,2)=c(1)*v(1208)
dcsdc2(1,2,3)=v(1283)
dcsdc2(1,2,4)=v(1284)
dcsdc2(1,2,5)=v(1285)
dcsdc2(1,2,6)=v(1286)
dcsdc2(1,2,7)=v(1287)
dcsdc2(1,2,8)=v(1288)
dcsdc2(1,2,9)=v(1289)
dcsdc2(1,2,10)=v(1290)
dcsdc2(1,2,11)=v(1291)
dcsdc2(1,2,12)=v(1292)
dcsdc2(1,3,1)=v(1212)
dcsdc2(1,3,2)=v(1283)
dcsdc2(1,3,3)=c(1)*v(1206)
dcsdc2(1,3,4)=v(1349)
dcsdc2(1,3,5)=v(1350)
dcsdc2(1,3,6)=v(1351)
dcsdc2(1,3,7)=v(1352)
dcsdc2(1,3,8)=v(1353)
dcsdc2(1,3,9)=v(1354)
dcsdc2(1,3,10)=v(1355)
dcsdc2(1,3,11)=v(1356)
dcsdc2(1,3,12)=v(1357)
dcsdc2(1,4,1)=v(1213)
dcsdc2(1,4,2)=v(1284)
dcsdc2(1,4,3)=v(1349)
dcsdc2(1,4,4)=c(1)*v(1203)
dcsdc2(1,4,5)=v(1409)
dcsdc2(1,4,6)=v(1410)
dcsdc2(1,4,7)=v(1411)
dcsdc2(1,4,8)=v(1412)
dcsdc2(1,4,9)=v(1413)
dcsdc2(1,4,10)=v(1414)
dcsdc2(1,4,11)=v(1415)
dcsdc2(1,4,12)=v(1416)
dcsdc2(1,5,1)=v(1214)
dcsdc2(1,5,2)=v(1285)
dcsdc2(1,5,3)=v(1350)
dcsdc2(1,5,4)=v(1409)
dcsdc2(1,5,5)=c(1)*v(1198)
dcsdc2(1,5,6)=v(1463)
dcsdc2(1,5,7)=v(1464)
dcsdc2(1,5,8)=v(1465)
dcsdc2(1,5,9)=v(1466)
dcsdc2(1,5,10)=v(1467)
dcsdc2(1,5,11)=v(1468)
dcsdc2(1,5,12)=v(1469)
dcsdc2(1,6,1)=v(1215)
dcsdc2(1,6,2)=v(1286)
dcsdc2(1,6,3)=v(1351)
dcsdc2(1,6,4)=v(1410)
dcsdc2(1,6,5)=v(1463)
dcsdc2(1,6,6)=c(1)*v(1191)
dcsdc2(1,6,7)=v(1511)
dcsdc2(1,6,8)=v(1512)
dcsdc2(1,6,9)=v(1513)
dcsdc2(1,6,10)=v(1514)
dcsdc2(1,6,11)=v(1515)
dcsdc2(1,6,12)=v(1516)
dcsdc2(1,7,1)=v(1216)
dcsdc2(1,7,2)=v(1287)
dcsdc2(1,7,3)=v(1352)
dcsdc2(1,7,4)=v(1411)
dcsdc2(1,7,5)=v(1464)
dcsdc2(1,7,6)=v(1511)
dcsdc2(1,7,7)=c(1)*v(1170)
dcsdc2(1,7,8)=v(1553)
dcsdc2(1,7,9)=v(1554)
dcsdc2(1,7,10)=v(1555)
dcsdc2(1,7,11)=v(1556)
dcsdc2(1,7,12)=v(1557)
dcsdc2(1,8,1)=v(1217)
dcsdc2(1,8,2)=v(1288)
dcsdc2(1,8,3)=v(1353)
dcsdc2(1,8,4)=v(1412)
dcsdc2(1,8,5)=v(1465)
dcsdc2(1,8,6)=v(1512)
dcsdc2(1,8,7)=v(1553)
dcsdc2(1,8,8)=c(1)*v(1163)
dcsdc2(1,8,9)=v(1589)
dcsdc2(1,8,10)=v(1590)
dcsdc2(1,8,11)=v(1591)
dcsdc2(1,8,12)=v(1592)
dcsdc2(1,9,1)=v(1218)
dcsdc2(1,9,2)=v(1289)
dcsdc2(1,9,3)=v(1354)
dcsdc2(1,9,4)=v(1413)
dcsdc2(1,9,5)=v(1466)
dcsdc2(1,9,6)=v(1513)
dcsdc2(1,9,7)=v(1554)
dcsdc2(1,9,8)=v(1589)
dcsdc2(1,9,9)=c(1)*v(1155)
dcsdc2(1,9,10)=v(1619)
dcsdc2(1,9,11)=v(1620)
dcsdc2(1,9,12)=v(1621)
dcsdc2(1,10,1)=v(1219)
dcsdc2(1,10,2)=v(1290)
dcsdc2(1,10,3)=v(1355)
dcsdc2(1,10,4)=v(1414)
dcsdc2(1,10,5)=v(1467)
dcsdc2(1,10,6)=v(1514)
dcsdc2(1,10,7)=v(1555)
dcsdc2(1,10,8)=v(1590)
dcsdc2(1,10,9)=v(1619)
dcsdc2(1,10,10)=c(1)*v(1146)
dcsdc2(1,10,11)=v(1643)
dcsdc2(1,10,12)=v(1644)
dcsdc2(1,11,1)=v(1220)
dcsdc2(1,11,2)=v(1291)
dcsdc2(1,11,3)=v(1356)
dcsdc2(1,11,4)=v(1415)
dcsdc2(1,11,5)=v(1468)
dcsdc2(1,11,6)=v(1515)
dcsdc2(1,11,7)=v(1556)
dcsdc2(1,11,8)=v(1591)
dcsdc2(1,11,9)=v(1620)
dcsdc2(1,11,10)=v(1643)
dcsdc2(1,11,11)=c(1)*v(1135)
dcsdc2(1,11,12)=v(1661)
dcsdc2(1,12,1)=v(1221)
dcsdc2(1,12,2)=v(1292)
dcsdc2(1,12,3)=v(1357)
dcsdc2(1,12,4)=v(1416)
dcsdc2(1,12,5)=v(1469)
dcsdc2(1,12,6)=v(1516)
dcsdc2(1,12,7)=v(1557)
dcsdc2(1,12,8)=v(1592)
dcsdc2(1,12,9)=v(1621)
dcsdc2(1,12,10)=v(1644)
dcsdc2(1,12,11)=v(1661)
dcsdc2(1,12,12)=c(1)*v(1122)
dcsdc2(2,1,1)=c(2)*v(1209)
dcsdc2(2,1,2)=v(1223)
dcsdc2(2,1,3)=v(1224)
dcsdc2(2,1,4)=v(1225)
dcsdc2(2,1,5)=v(1226)
dcsdc2(2,1,6)=v(1227)
dcsdc2(2,1,7)=v(1228)
dcsdc2(2,1,8)=v(1229)
dcsdc2(2,1,9)=v(1230)
dcsdc2(2,1,10)=v(1231)
dcsdc2(2,1,11)=v(1232)
dcsdc2(2,1,12)=v(1233)
dcsdc2(2,2,1)=v(1223)
dcsdc2(2,2,2)=c(2)*v(1208)+2d0*v(968)
dcsdc2(2,2,3)=v(1294)
dcsdc2(2,2,4)=v(1295)
dcsdc2(2,2,5)=v(1296)
dcsdc2(2,2,6)=v(1297)
dcsdc2(2,2,7)=v(1298)
dcsdc2(2,2,8)=v(1299)
dcsdc2(2,2,9)=v(1300)
dcsdc2(2,2,10)=v(1301)
dcsdc2(2,2,11)=v(1302)
dcsdc2(2,2,12)=v(1303)
dcsdc2(2,3,1)=v(1224)
dcsdc2(2,3,2)=v(1294)
dcsdc2(2,3,3)=c(2)*v(1206)
dcsdc2(2,3,4)=v(1359)
dcsdc2(2,3,5)=v(1360)
dcsdc2(2,3,6)=v(1361)
dcsdc2(2,3,7)=v(1362)
dcsdc2(2,3,8)=v(1363)
dcsdc2(2,3,9)=v(1364)
dcsdc2(2,3,10)=v(1365)
dcsdc2(2,3,11)=v(1366)
dcsdc2(2,3,12)=v(1367)
dcsdc2(2,4,1)=v(1225)
dcsdc2(2,4,2)=v(1295)
dcsdc2(2,4,3)=v(1359)
dcsdc2(2,4,4)=c(2)*v(1203)
dcsdc2(2,4,5)=v(1418)
dcsdc2(2,4,6)=v(1419)
dcsdc2(2,4,7)=v(1420)
dcsdc2(2,4,8)=v(1421)
dcsdc2(2,4,9)=v(1422)
dcsdc2(2,4,10)=v(1423)
dcsdc2(2,4,11)=v(1424)
dcsdc2(2,4,12)=v(1425)
dcsdc2(2,5,1)=v(1226)
dcsdc2(2,5,2)=v(1296)
dcsdc2(2,5,3)=v(1360)
dcsdc2(2,5,4)=v(1418)
dcsdc2(2,5,5)=c(2)*v(1198)
dcsdc2(2,5,6)=v(1471)
dcsdc2(2,5,7)=v(1472)
dcsdc2(2,5,8)=v(1473)
dcsdc2(2,5,9)=v(1474)
dcsdc2(2,5,10)=v(1475)
dcsdc2(2,5,11)=v(1476)
dcsdc2(2,5,12)=v(1477)
dcsdc2(2,6,1)=v(1227)
dcsdc2(2,6,2)=v(1297)
dcsdc2(2,6,3)=v(1361)
dcsdc2(2,6,4)=v(1419)
dcsdc2(2,6,5)=v(1471)
dcsdc2(2,6,6)=c(2)*v(1191)
dcsdc2(2,6,7)=v(1518)
dcsdc2(2,6,8)=v(1519)
dcsdc2(2,6,9)=v(1520)
dcsdc2(2,6,10)=v(1521)
dcsdc2(2,6,11)=v(1522)
dcsdc2(2,6,12)=v(1523)
dcsdc2(2,7,1)=v(1228)
dcsdc2(2,7,2)=v(1298)
dcsdc2(2,7,3)=v(1362)
dcsdc2(2,7,4)=v(1420)
dcsdc2(2,7,5)=v(1472)
dcsdc2(2,7,6)=v(1518)
dcsdc2(2,7,7)=c(2)*v(1170)
dcsdc2(2,7,8)=v(1559)
dcsdc2(2,7,9)=v(1560)
dcsdc2(2,7,10)=v(1561)
dcsdc2(2,7,11)=v(1562)
dcsdc2(2,7,12)=v(1563)
dcsdc2(2,8,1)=v(1229)
dcsdc2(2,8,2)=v(1299)
dcsdc2(2,8,3)=v(1363)
dcsdc2(2,8,4)=v(1421)
dcsdc2(2,8,5)=v(1473)
dcsdc2(2,8,6)=v(1519)
dcsdc2(2,8,7)=v(1559)
dcsdc2(2,8,8)=c(2)*v(1163)
dcsdc2(2,8,9)=v(1594)
dcsdc2(2,8,10)=v(1595)
dcsdc2(2,8,11)=v(1596)
dcsdc2(2,8,12)=v(1597)
dcsdc2(2,9,1)=v(1230)
dcsdc2(2,9,2)=v(1300)
dcsdc2(2,9,3)=v(1364)
dcsdc2(2,9,4)=v(1422)
dcsdc2(2,9,5)=v(1474)
dcsdc2(2,9,6)=v(1520)
dcsdc2(2,9,7)=v(1560)
dcsdc2(2,9,8)=v(1594)
dcsdc2(2,9,9)=c(2)*v(1155)
dcsdc2(2,9,10)=v(1623)
dcsdc2(2,9,11)=v(1624)
dcsdc2(2,9,12)=v(1625)
dcsdc2(2,10,1)=v(1231)
dcsdc2(2,10,2)=v(1301)
dcsdc2(2,10,3)=v(1365)
dcsdc2(2,10,4)=v(1423)
dcsdc2(2,10,5)=v(1475)
dcsdc2(2,10,6)=v(1521)
dcsdc2(2,10,7)=v(1561)
dcsdc2(2,10,8)=v(1595)
dcsdc2(2,10,9)=v(1623)
dcsdc2(2,10,10)=c(2)*v(1146)
dcsdc2(2,10,11)=v(1646)
dcsdc2(2,10,12)=v(1647)
dcsdc2(2,11,1)=v(1232)
dcsdc2(2,11,2)=v(1302)
dcsdc2(2,11,3)=v(1366)
dcsdc2(2,11,4)=v(1424)
dcsdc2(2,11,5)=v(1476)
dcsdc2(2,11,6)=v(1522)
dcsdc2(2,11,7)=v(1562)
dcsdc2(2,11,8)=v(1596)
dcsdc2(2,11,9)=v(1624)
dcsdc2(2,11,10)=v(1646)
dcsdc2(2,11,11)=c(2)*v(1135)
dcsdc2(2,11,12)=v(1663)
dcsdc2(2,12,1)=v(1233)
dcsdc2(2,12,2)=v(1303)
dcsdc2(2,12,3)=v(1367)
dcsdc2(2,12,4)=v(1425)
dcsdc2(2,12,5)=v(1477)
dcsdc2(2,12,6)=v(1523)
dcsdc2(2,12,7)=v(1563)
dcsdc2(2,12,8)=v(1597)
dcsdc2(2,12,9)=v(1625)
dcsdc2(2,12,10)=v(1647)
dcsdc2(2,12,11)=v(1663)
dcsdc2(2,12,12)=c(2)*v(1122)
dcsdc2(3,1,1)=c(3)*v(1209)
dcsdc2(3,1,2)=v(1235)
dcsdc2(3,1,3)=v(1236)
dcsdc2(3,1,4)=v(1237)
dcsdc2(3,1,5)=v(1238)
dcsdc2(3,1,6)=v(1239)
dcsdc2(3,1,7)=v(1240)
dcsdc2(3,1,8)=v(1241)
dcsdc2(3,1,9)=v(1242)
dcsdc2(3,1,10)=v(1243)
dcsdc2(3,1,11)=v(1244)
dcsdc2(3,1,12)=v(1245)
dcsdc2(3,2,1)=v(1235)
dcsdc2(3,2,2)=c(3)*v(1208)
dcsdc2(3,2,3)=v(1305)
dcsdc2(3,2,4)=v(1306)
dcsdc2(3,2,5)=v(1307)
dcsdc2(3,2,6)=v(1308)
dcsdc2(3,2,7)=v(1309)
dcsdc2(3,2,8)=v(1310)
dcsdc2(3,2,9)=v(1311)
dcsdc2(3,2,10)=v(1312)
dcsdc2(3,2,11)=v(1313)
dcsdc2(3,2,12)=v(1314)
dcsdc2(3,3,1)=v(1236)
dcsdc2(3,3,2)=v(1305)
dcsdc2(3,3,3)=c(3)*v(1206)+2d0*v(969)
dcsdc2(3,3,4)=v(1369)
dcsdc2(3,3,5)=v(1370)
dcsdc2(3,3,6)=v(1371)
dcsdc2(3,3,7)=v(1372)
dcsdc2(3,3,8)=v(1373)
dcsdc2(3,3,9)=v(1374)
dcsdc2(3,3,10)=v(1375)
dcsdc2(3,3,11)=v(1376)
dcsdc2(3,3,12)=v(1377)
dcsdc2(3,4,1)=v(1237)
dcsdc2(3,4,2)=v(1306)
dcsdc2(3,4,3)=v(1369)
dcsdc2(3,4,4)=c(3)*v(1203)
dcsdc2(3,4,5)=v(1427)
dcsdc2(3,4,6)=v(1428)
dcsdc2(3,4,7)=v(1429)
dcsdc2(3,4,8)=v(1430)
dcsdc2(3,4,9)=v(1431)
dcsdc2(3,4,10)=v(1432)
dcsdc2(3,4,11)=v(1433)
dcsdc2(3,4,12)=v(1434)
dcsdc2(3,5,1)=v(1238)
dcsdc2(3,5,2)=v(1307)
dcsdc2(3,5,3)=v(1370)
dcsdc2(3,5,4)=v(1427)
dcsdc2(3,5,5)=c(3)*v(1198)
dcsdc2(3,5,6)=v(1479)
dcsdc2(3,5,7)=v(1480)
dcsdc2(3,5,8)=v(1481)
dcsdc2(3,5,9)=v(1482)
dcsdc2(3,5,10)=v(1483)
dcsdc2(3,5,11)=v(1484)
dcsdc2(3,5,12)=v(1485)
dcsdc2(3,6,1)=v(1239)
dcsdc2(3,6,2)=v(1308)
dcsdc2(3,6,3)=v(1371)
dcsdc2(3,6,4)=v(1428)
dcsdc2(3,6,5)=v(1479)
dcsdc2(3,6,6)=c(3)*v(1191)
dcsdc2(3,6,7)=v(1525)
dcsdc2(3,6,8)=v(1526)
dcsdc2(3,6,9)=v(1527)
dcsdc2(3,6,10)=v(1528)
dcsdc2(3,6,11)=v(1529)
dcsdc2(3,6,12)=v(1530)
dcsdc2(3,7,1)=v(1240)
dcsdc2(3,7,2)=v(1309)
dcsdc2(3,7,3)=v(1372)
dcsdc2(3,7,4)=v(1429)
dcsdc2(3,7,5)=v(1480)
dcsdc2(3,7,6)=v(1525)
dcsdc2(3,7,7)=c(3)*v(1170)
dcsdc2(3,7,8)=v(1565)
dcsdc2(3,7,9)=v(1566)
dcsdc2(3,7,10)=v(1567)
dcsdc2(3,7,11)=v(1568)
dcsdc2(3,7,12)=v(1569)
dcsdc2(3,8,1)=v(1241)
dcsdc2(3,8,2)=v(1310)
dcsdc2(3,8,3)=v(1373)
dcsdc2(3,8,4)=v(1430)
dcsdc2(3,8,5)=v(1481)
dcsdc2(3,8,6)=v(1526)
dcsdc2(3,8,7)=v(1565)
dcsdc2(3,8,8)=c(3)*v(1163)
dcsdc2(3,8,9)=v(1599)
dcsdc2(3,8,10)=v(1600)
dcsdc2(3,8,11)=v(1601)
dcsdc2(3,8,12)=v(1602)
dcsdc2(3,9,1)=v(1242)
dcsdc2(3,9,2)=v(1311)
dcsdc2(3,9,3)=v(1374)
dcsdc2(3,9,4)=v(1431)
dcsdc2(3,9,5)=v(1482)
dcsdc2(3,9,6)=v(1527)
dcsdc2(3,9,7)=v(1566)
dcsdc2(3,9,8)=v(1599)
dcsdc2(3,9,9)=c(3)*v(1155)
dcsdc2(3,9,10)=v(1627)
dcsdc2(3,9,11)=v(1628)
dcsdc2(3,9,12)=v(1629)
dcsdc2(3,10,1)=v(1243)
dcsdc2(3,10,2)=v(1312)
dcsdc2(3,10,3)=v(1375)
dcsdc2(3,10,4)=v(1432)
dcsdc2(3,10,5)=v(1483)
dcsdc2(3,10,6)=v(1528)
dcsdc2(3,10,7)=v(1567)
dcsdc2(3,10,8)=v(1600)
dcsdc2(3,10,9)=v(1627)
dcsdc2(3,10,10)=c(3)*v(1146)
dcsdc2(3,10,11)=v(1649)
dcsdc2(3,10,12)=v(1650)
dcsdc2(3,11,1)=v(1244)
dcsdc2(3,11,2)=v(1313)
dcsdc2(3,11,3)=v(1376)
dcsdc2(3,11,4)=v(1433)
dcsdc2(3,11,5)=v(1484)
dcsdc2(3,11,6)=v(1529)
dcsdc2(3,11,7)=v(1568)
dcsdc2(3,11,8)=v(1601)
dcsdc2(3,11,9)=v(1628)
dcsdc2(3,11,10)=v(1649)
dcsdc2(3,11,11)=c(3)*v(1135)
dcsdc2(3,11,12)=v(1665)
dcsdc2(3,12,1)=v(1245)
dcsdc2(3,12,2)=v(1314)
dcsdc2(3,12,3)=v(1377)
dcsdc2(3,12,4)=v(1434)
dcsdc2(3,12,5)=v(1485)
dcsdc2(3,12,6)=v(1530)
dcsdc2(3,12,7)=v(1569)
dcsdc2(3,12,8)=v(1602)
dcsdc2(3,12,9)=v(1629)
dcsdc2(3,12,10)=v(1650)
dcsdc2(3,12,11)=v(1665)
dcsdc2(3,12,12)=c(3)*v(1122)
dcsdc2(4,1,1)=c(4)*v(1209)
dcsdc2(4,1,2)=v(1247)
dcsdc2(4,1,3)=v(1248)
dcsdc2(4,1,4)=v(1249)
dcsdc2(4,1,5)=v(1250)
dcsdc2(4,1,6)=v(1251)
dcsdc2(4,1,7)=v(1252)
dcsdc2(4,1,8)=v(1253)
dcsdc2(4,1,9)=v(1254)
dcsdc2(4,1,10)=v(1255)
dcsdc2(4,1,11)=v(1256)
dcsdc2(4,1,12)=v(1257)
dcsdc2(4,2,1)=v(1247)
dcsdc2(4,2,2)=c(4)*v(1208)
dcsdc2(4,2,3)=v(1316)
dcsdc2(4,2,4)=v(1317)
dcsdc2(4,2,5)=v(1318)
dcsdc2(4,2,6)=v(1319)
dcsdc2(4,2,7)=v(1320)
dcsdc2(4,2,8)=v(1321)
dcsdc2(4,2,9)=v(1322)
dcsdc2(4,2,10)=v(1323)
dcsdc2(4,2,11)=v(1324)
dcsdc2(4,2,12)=v(1325)
dcsdc2(4,3,1)=v(1248)
dcsdc2(4,3,2)=v(1316)
dcsdc2(4,3,3)=c(4)*v(1206)
dcsdc2(4,3,4)=v(1379)
dcsdc2(4,3,5)=v(1380)
dcsdc2(4,3,6)=v(1381)
dcsdc2(4,3,7)=v(1382)
dcsdc2(4,3,8)=v(1383)
dcsdc2(4,3,9)=v(1384)
dcsdc2(4,3,10)=v(1385)
dcsdc2(4,3,11)=v(1386)
dcsdc2(4,3,12)=v(1387)
dcsdc2(4,4,1)=v(1249)
dcsdc2(4,4,2)=v(1317)
dcsdc2(4,4,3)=v(1379)
dcsdc2(4,4,4)=c(4)*v(1203)+2d0*v(970)
dcsdc2(4,4,5)=v(1436)
dcsdc2(4,4,6)=v(1437)
dcsdc2(4,4,7)=v(1438)
dcsdc2(4,4,8)=v(1439)
dcsdc2(4,4,9)=v(1440)
dcsdc2(4,4,10)=v(1441)
dcsdc2(4,4,11)=v(1442)
dcsdc2(4,4,12)=v(1443)
dcsdc2(4,5,1)=v(1250)
dcsdc2(4,5,2)=v(1318)
dcsdc2(4,5,3)=v(1380)
dcsdc2(4,5,4)=v(1436)
dcsdc2(4,5,5)=c(4)*v(1198)
dcsdc2(4,5,6)=v(1487)
dcsdc2(4,5,7)=v(1488)
dcsdc2(4,5,8)=v(1489)
dcsdc2(4,5,9)=v(1490)
dcsdc2(4,5,10)=v(1491)
dcsdc2(4,5,11)=v(1492)
dcsdc2(4,5,12)=v(1493)
dcsdc2(4,6,1)=v(1251)
dcsdc2(4,6,2)=v(1319)
dcsdc2(4,6,3)=v(1381)
dcsdc2(4,6,4)=v(1437)
dcsdc2(4,6,5)=v(1487)
dcsdc2(4,6,6)=c(4)*v(1191)
dcsdc2(4,6,7)=v(1532)
dcsdc2(4,6,8)=v(1533)
dcsdc2(4,6,9)=v(1534)
dcsdc2(4,6,10)=v(1535)
dcsdc2(4,6,11)=v(1536)
dcsdc2(4,6,12)=v(1537)
dcsdc2(4,7,1)=v(1252)
dcsdc2(4,7,2)=v(1320)
dcsdc2(4,7,3)=v(1382)
dcsdc2(4,7,4)=v(1438)
dcsdc2(4,7,5)=v(1488)
dcsdc2(4,7,6)=v(1532)
dcsdc2(4,7,7)=c(4)*v(1170)
dcsdc2(4,7,8)=v(1571)
dcsdc2(4,7,9)=v(1572)
dcsdc2(4,7,10)=v(1573)
dcsdc2(4,7,11)=v(1574)
dcsdc2(4,7,12)=v(1575)
dcsdc2(4,8,1)=v(1253)
dcsdc2(4,8,2)=v(1321)
dcsdc2(4,8,3)=v(1383)
dcsdc2(4,8,4)=v(1439)
dcsdc2(4,8,5)=v(1489)
dcsdc2(4,8,6)=v(1533)
dcsdc2(4,8,7)=v(1571)
dcsdc2(4,8,8)=c(4)*v(1163)
dcsdc2(4,8,9)=v(1604)
dcsdc2(4,8,10)=v(1605)
dcsdc2(4,8,11)=v(1606)
dcsdc2(4,8,12)=v(1607)
dcsdc2(4,9,1)=v(1254)
dcsdc2(4,9,2)=v(1322)
dcsdc2(4,9,3)=v(1384)
dcsdc2(4,9,4)=v(1440)
dcsdc2(4,9,5)=v(1490)
dcsdc2(4,9,6)=v(1534)
dcsdc2(4,9,7)=v(1572)
dcsdc2(4,9,8)=v(1604)
dcsdc2(4,9,9)=c(4)*v(1155)
dcsdc2(4,9,10)=v(1631)
dcsdc2(4,9,11)=v(1632)
dcsdc2(4,9,12)=v(1633)
dcsdc2(4,10,1)=v(1255)
dcsdc2(4,10,2)=v(1323)
dcsdc2(4,10,3)=v(1385)
dcsdc2(4,10,4)=v(1441)
dcsdc2(4,10,5)=v(1491)
dcsdc2(4,10,6)=v(1535)
dcsdc2(4,10,7)=v(1573)
dcsdc2(4,10,8)=v(1605)
dcsdc2(4,10,9)=v(1631)
dcsdc2(4,10,10)=c(4)*v(1146)
dcsdc2(4,10,11)=v(1652)
dcsdc2(4,10,12)=v(1653)
dcsdc2(4,11,1)=v(1256)
dcsdc2(4,11,2)=v(1324)
dcsdc2(4,11,3)=v(1386)
dcsdc2(4,11,4)=v(1442)
dcsdc2(4,11,5)=v(1492)
dcsdc2(4,11,6)=v(1536)
dcsdc2(4,11,7)=v(1574)
dcsdc2(4,11,8)=v(1606)
dcsdc2(4,11,9)=v(1632)
dcsdc2(4,11,10)=v(1652)
dcsdc2(4,11,11)=c(4)*v(1135)
dcsdc2(4,11,12)=v(1667)
dcsdc2(4,12,1)=v(1257)
dcsdc2(4,12,2)=v(1325)
dcsdc2(4,12,3)=v(1387)
dcsdc2(4,12,4)=v(1443)
dcsdc2(4,12,5)=v(1493)
dcsdc2(4,12,6)=v(1537)
dcsdc2(4,12,7)=v(1575)
dcsdc2(4,12,8)=v(1607)
dcsdc2(4,12,9)=v(1633)
dcsdc2(4,12,10)=v(1653)
dcsdc2(4,12,11)=v(1667)
dcsdc2(4,12,12)=c(4)*v(1122)
dcsdc2(5,1,1)=c(5)*v(1209)
dcsdc2(5,1,2)=v(1259)
dcsdc2(5,1,3)=v(1260)
dcsdc2(5,1,4)=v(1261)
dcsdc2(5,1,5)=v(1262)
dcsdc2(5,1,6)=v(1263)
dcsdc2(5,1,7)=v(1264)
dcsdc2(5,1,8)=v(1265)
dcsdc2(5,1,9)=v(1266)
dcsdc2(5,1,10)=v(1267)
dcsdc2(5,1,11)=v(1268)
dcsdc2(5,1,12)=v(1269)
dcsdc2(5,2,1)=v(1259)
dcsdc2(5,2,2)=c(5)*v(1208)
dcsdc2(5,2,3)=v(1327)
dcsdc2(5,2,4)=v(1328)
dcsdc2(5,2,5)=v(1329)
dcsdc2(5,2,6)=v(1330)
dcsdc2(5,2,7)=v(1331)
dcsdc2(5,2,8)=v(1332)
dcsdc2(5,2,9)=v(1333)
dcsdc2(5,2,10)=v(1334)
dcsdc2(5,2,11)=v(1335)
dcsdc2(5,2,12)=v(1336)
dcsdc2(5,3,1)=v(1260)
dcsdc2(5,3,2)=v(1327)
dcsdc2(5,3,3)=c(5)*v(1206)
dcsdc2(5,3,4)=v(1389)
dcsdc2(5,3,5)=v(1390)
dcsdc2(5,3,6)=v(1391)
dcsdc2(5,3,7)=v(1392)
dcsdc2(5,3,8)=v(1393)
dcsdc2(5,3,9)=v(1394)
dcsdc2(5,3,10)=v(1395)
dcsdc2(5,3,11)=v(1396)
dcsdc2(5,3,12)=v(1397)
dcsdc2(5,4,1)=v(1261)
dcsdc2(5,4,2)=v(1328)
dcsdc2(5,4,3)=v(1389)
dcsdc2(5,4,4)=c(5)*v(1203)
dcsdc2(5,4,5)=v(1445)
dcsdc2(5,4,6)=v(1446)
dcsdc2(5,4,7)=v(1447)
dcsdc2(5,4,8)=v(1448)
dcsdc2(5,4,9)=v(1449)
dcsdc2(5,4,10)=v(1450)
dcsdc2(5,4,11)=v(1451)
dcsdc2(5,4,12)=v(1452)
dcsdc2(5,5,1)=v(1262)
dcsdc2(5,5,2)=v(1329)
dcsdc2(5,5,3)=v(1390)
dcsdc2(5,5,4)=v(1445)
dcsdc2(5,5,5)=c(5)*v(1198)+2d0*v(971)
dcsdc2(5,5,6)=v(1495)
dcsdc2(5,5,7)=v(1496)
dcsdc2(5,5,8)=v(1497)
dcsdc2(5,5,9)=v(1498)
dcsdc2(5,5,10)=v(1499)
dcsdc2(5,5,11)=v(1500)
dcsdc2(5,5,12)=v(1501)
dcsdc2(5,6,1)=v(1263)
dcsdc2(5,6,2)=v(1330)
dcsdc2(5,6,3)=v(1391)
dcsdc2(5,6,4)=v(1446)
dcsdc2(5,6,5)=v(1495)
dcsdc2(5,6,6)=c(5)*v(1191)
dcsdc2(5,6,7)=v(1539)
dcsdc2(5,6,8)=v(1540)
dcsdc2(5,6,9)=v(1541)
dcsdc2(5,6,10)=v(1542)
dcsdc2(5,6,11)=v(1543)
dcsdc2(5,6,12)=v(1544)
dcsdc2(5,7,1)=v(1264)
dcsdc2(5,7,2)=v(1331)
dcsdc2(5,7,3)=v(1392)
dcsdc2(5,7,4)=v(1447)
dcsdc2(5,7,5)=v(1496)
dcsdc2(5,7,6)=v(1539)
dcsdc2(5,7,7)=c(5)*v(1170)
dcsdc2(5,7,8)=v(1577)
dcsdc2(5,7,9)=v(1578)
dcsdc2(5,7,10)=v(1579)
dcsdc2(5,7,11)=v(1580)
dcsdc2(5,7,12)=v(1581)
dcsdc2(5,8,1)=v(1265)
dcsdc2(5,8,2)=v(1332)
dcsdc2(5,8,3)=v(1393)
dcsdc2(5,8,4)=v(1448)
dcsdc2(5,8,5)=v(1497)
dcsdc2(5,8,6)=v(1540)
dcsdc2(5,8,7)=v(1577)
dcsdc2(5,8,8)=c(5)*v(1163)
dcsdc2(5,8,9)=v(1609)
dcsdc2(5,8,10)=v(1610)
dcsdc2(5,8,11)=v(1611)
dcsdc2(5,8,12)=v(1612)
dcsdc2(5,9,1)=v(1266)
dcsdc2(5,9,2)=v(1333)
dcsdc2(5,9,3)=v(1394)
dcsdc2(5,9,4)=v(1449)
dcsdc2(5,9,5)=v(1498)
dcsdc2(5,9,6)=v(1541)
dcsdc2(5,9,7)=v(1578)
dcsdc2(5,9,8)=v(1609)
dcsdc2(5,9,9)=c(5)*v(1155)
dcsdc2(5,9,10)=v(1635)
dcsdc2(5,9,11)=v(1636)
dcsdc2(5,9,12)=v(1637)
dcsdc2(5,10,1)=v(1267)
dcsdc2(5,10,2)=v(1334)
dcsdc2(5,10,3)=v(1395)
dcsdc2(5,10,4)=v(1450)
dcsdc2(5,10,5)=v(1499)
dcsdc2(5,10,6)=v(1542)
dcsdc2(5,10,7)=v(1579)
dcsdc2(5,10,8)=v(1610)
dcsdc2(5,10,9)=v(1635)
dcsdc2(5,10,10)=c(5)*v(1146)
dcsdc2(5,10,11)=v(1655)
dcsdc2(5,10,12)=v(1656)
dcsdc2(5,11,1)=v(1268)
dcsdc2(5,11,2)=v(1335)
dcsdc2(5,11,3)=v(1396)
dcsdc2(5,11,4)=v(1451)
dcsdc2(5,11,5)=v(1500)
dcsdc2(5,11,6)=v(1543)
dcsdc2(5,11,7)=v(1580)
dcsdc2(5,11,8)=v(1611)
dcsdc2(5,11,9)=v(1636)
dcsdc2(5,11,10)=v(1655)
dcsdc2(5,11,11)=c(5)*v(1135)
dcsdc2(5,11,12)=v(1669)
dcsdc2(5,12,1)=v(1269)
dcsdc2(5,12,2)=v(1336)
dcsdc2(5,12,3)=v(1397)
dcsdc2(5,12,4)=v(1452)
dcsdc2(5,12,5)=v(1501)
dcsdc2(5,12,6)=v(1544)
dcsdc2(5,12,7)=v(1581)
dcsdc2(5,12,8)=v(1612)
dcsdc2(5,12,9)=v(1637)
dcsdc2(5,12,10)=v(1656)
dcsdc2(5,12,11)=v(1669)
dcsdc2(5,12,12)=c(5)*v(1122)
dcsdc2(6,1,1)=c(6)*v(1209)
dcsdc2(6,1,2)=v(1271)
dcsdc2(6,1,3)=v(1272)
dcsdc2(6,1,4)=v(1273)
dcsdc2(6,1,5)=v(1274)
dcsdc2(6,1,6)=v(1275)
dcsdc2(6,1,7)=v(1276)
dcsdc2(6,1,8)=v(1277)
dcsdc2(6,1,9)=v(1278)
dcsdc2(6,1,10)=v(1279)
dcsdc2(6,1,11)=v(1280)
dcsdc2(6,1,12)=v(1281)
dcsdc2(6,2,1)=v(1271)
dcsdc2(6,2,2)=c(6)*v(1208)
dcsdc2(6,2,3)=v(1338)
dcsdc2(6,2,4)=v(1339)
dcsdc2(6,2,5)=v(1340)
dcsdc2(6,2,6)=v(1341)
dcsdc2(6,2,7)=v(1342)
dcsdc2(6,2,8)=v(1343)
dcsdc2(6,2,9)=v(1344)
dcsdc2(6,2,10)=v(1345)
dcsdc2(6,2,11)=v(1346)
dcsdc2(6,2,12)=v(1347)
dcsdc2(6,3,1)=v(1272)
dcsdc2(6,3,2)=v(1338)
dcsdc2(6,3,3)=c(6)*v(1206)
dcsdc2(6,3,4)=v(1399)
dcsdc2(6,3,5)=v(1400)
dcsdc2(6,3,6)=v(1401)
dcsdc2(6,3,7)=v(1402)
dcsdc2(6,3,8)=v(1403)
dcsdc2(6,3,9)=v(1404)
dcsdc2(6,3,10)=v(1405)
dcsdc2(6,3,11)=v(1406)
dcsdc2(6,3,12)=v(1407)
dcsdc2(6,4,1)=v(1273)
dcsdc2(6,4,2)=v(1339)
dcsdc2(6,4,3)=v(1399)
dcsdc2(6,4,4)=c(6)*v(1203)
dcsdc2(6,4,5)=v(1454)
dcsdc2(6,4,6)=v(1455)
dcsdc2(6,4,7)=v(1456)
dcsdc2(6,4,8)=v(1457)
dcsdc2(6,4,9)=v(1458)
dcsdc2(6,4,10)=v(1459)
dcsdc2(6,4,11)=v(1460)
dcsdc2(6,4,12)=v(1461)
dcsdc2(6,5,1)=v(1274)
dcsdc2(6,5,2)=v(1340)
dcsdc2(6,5,3)=v(1400)
dcsdc2(6,5,4)=v(1454)
dcsdc2(6,5,5)=c(6)*v(1198)
dcsdc2(6,5,6)=v(1503)
dcsdc2(6,5,7)=v(1504)
dcsdc2(6,5,8)=v(1505)
dcsdc2(6,5,9)=v(1506)
dcsdc2(6,5,10)=v(1507)
dcsdc2(6,5,11)=v(1508)
dcsdc2(6,5,12)=v(1509)
dcsdc2(6,6,1)=v(1275)
dcsdc2(6,6,2)=v(1341)
dcsdc2(6,6,3)=v(1401)
dcsdc2(6,6,4)=v(1455)
dcsdc2(6,6,5)=v(1503)
dcsdc2(6,6,6)=c(6)*v(1191)+2d0*v(972)
dcsdc2(6,6,7)=v(1546)
dcsdc2(6,6,8)=v(1547)
dcsdc2(6,6,9)=v(1548)
dcsdc2(6,6,10)=v(1549)
dcsdc2(6,6,11)=v(1550)
dcsdc2(6,6,12)=v(1551)
dcsdc2(6,7,1)=v(1276)
dcsdc2(6,7,2)=v(1342)
dcsdc2(6,7,3)=v(1402)
dcsdc2(6,7,4)=v(1456)
dcsdc2(6,7,5)=v(1504)
dcsdc2(6,7,6)=v(1546)
dcsdc2(6,7,7)=c(6)*v(1170)
dcsdc2(6,7,8)=v(1583)
dcsdc2(6,7,9)=v(1584)
dcsdc2(6,7,10)=v(1585)
dcsdc2(6,7,11)=v(1586)
dcsdc2(6,7,12)=v(1587)
dcsdc2(6,8,1)=v(1277)
dcsdc2(6,8,2)=v(1343)
dcsdc2(6,8,3)=v(1403)
dcsdc2(6,8,4)=v(1457)
dcsdc2(6,8,5)=v(1505)
dcsdc2(6,8,6)=v(1547)
dcsdc2(6,8,7)=v(1583)
dcsdc2(6,8,8)=c(6)*v(1163)
dcsdc2(6,8,9)=v(1614)
dcsdc2(6,8,10)=v(1615)
dcsdc2(6,8,11)=v(1616)
dcsdc2(6,8,12)=v(1617)
dcsdc2(6,9,1)=v(1278)
dcsdc2(6,9,2)=v(1344)
dcsdc2(6,9,3)=v(1404)
dcsdc2(6,9,4)=v(1458)
dcsdc2(6,9,5)=v(1506)
dcsdc2(6,9,6)=v(1548)
dcsdc2(6,9,7)=v(1584)
dcsdc2(6,9,8)=v(1614)
dcsdc2(6,9,9)=c(6)*v(1155)
dcsdc2(6,9,10)=v(1639)
dcsdc2(6,9,11)=v(1640)
dcsdc2(6,9,12)=v(1641)
dcsdc2(6,10,1)=v(1279)
dcsdc2(6,10,2)=v(1345)
dcsdc2(6,10,3)=v(1405)
dcsdc2(6,10,4)=v(1459)
dcsdc2(6,10,5)=v(1507)
dcsdc2(6,10,6)=v(1549)
dcsdc2(6,10,7)=v(1585)
dcsdc2(6,10,8)=v(1615)
dcsdc2(6,10,9)=v(1639)
dcsdc2(6,10,10)=c(6)*v(1146)
dcsdc2(6,10,11)=v(1658)
dcsdc2(6,10,12)=v(1659)
dcsdc2(6,11,1)=v(1280)
dcsdc2(6,11,2)=v(1346)
dcsdc2(6,11,3)=v(1406)
dcsdc2(6,11,4)=v(1460)
dcsdc2(6,11,5)=v(1508)
dcsdc2(6,11,6)=v(1550)
dcsdc2(6,11,7)=v(1586)
dcsdc2(6,11,8)=v(1616)
dcsdc2(6,11,9)=v(1640)
dcsdc2(6,11,10)=v(1658)
dcsdc2(6,11,11)=c(6)*v(1135)
dcsdc2(6,11,12)=v(1671)
dcsdc2(6,12,1)=v(1281)
dcsdc2(6,12,2)=v(1347)
dcsdc2(6,12,3)=v(1407)
dcsdc2(6,12,4)=v(1461)
dcsdc2(6,12,5)=v(1509)
dcsdc2(6,12,6)=v(1551)
dcsdc2(6,12,7)=v(1587)
dcsdc2(6,12,8)=v(1617)
dcsdc2(6,12,9)=v(1641)
dcsdc2(6,12,10)=v(1659)
dcsdc2(6,12,11)=v(1671)
dcsdc2(6,12,12)=c(6)*v(1122)
END SUBROUTINE mksfbar
