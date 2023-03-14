//g++ sequential_intersect.cpp -o sequential_intersect -lgeos && ./sequential_intersect
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <geos/io/WKTReader.h>
#include <geos/io/WKTReader.inl>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Geometry.h>
#include <geos/geom/Envelope.h>
#include <geos/geom/Envelope.inl>
#include <geos/index/strtree/STRtree.h>
using namespace std;
using namespace geos::io;
using namespace geos::geom;
using namespace geos::index::strtree;
class candidate
{
public:
	int indexa; //�洢��ѡ�Ե�r�����ݼ�R�е��±� 
	int indexb; //�洢��ѡ�Ե�s�����ݼ�S�е��±�
	double centerx; //��ѡ�����ĵ��x���� 
	double centery; //��ѡ�����ĵ��y����
	double weight; //��ѡ�Ե�Ȩ�� 
};
vector<candidate> C; //�洢��ѡ�� 
vector<Geometry*> a; //�洢R����Ϣ 
vector<Geometry*> b; //�洢S����Ϣ 
clock_t astart;
void find_candidates()
{
	ifstream file1("WKT/STATE.csv");
	string wkt;
	STRtree treea; //R��R��
	GeometryFactory::Ptr factory = GeometryFactory::create();
	WKTReader reader(*factory); //��wktת��Ϊgeom
	int indexa = 0;
	astart = clock();
	while (getline(file1, wkt))
	{
		wkt = wkt.substr(1, wkt.rfind("\"") - 1); //csv�ļ�����wkt
		//wkt = wkt.substr(wkt.find("\t")+1, wkt.rfind("\t")-wkt.find("\t")-1); //�����ƻ�tsv�ļ�����wkt
		unique_ptr<Geometry> geom = reader.read(wkt); //����geom
		geom->setSRID(indexa); //����R��geom��SRID
		const Envelope* envelope = geom->getEnvelopeInternal(); //����MBR
		treea.insert(envelope, geom.get()); //����MBR
		a.push_back((Geometry*)geom.release());
		indexa++;
	}
	cout << (double)(clock() - astart) / CLOCKS_PER_SEC << endl;
	ifstream file2("WKT/PRIMARYROADS.csv");
	int indexb = 0;
	astart = clock();
	while (getline(file2, wkt))
	{
		wkt = wkt.substr(1, wkt.rfind("\"") - 1); //csv�ļ�����wkt
		//wkt = wkt.substr(wkt.find("\t")+1, wkt.rfind("\t")-wkt.find("\t")-1); //�����ƻ�tsv�ļ�����wkt
		unique_ptr<Geometry> geom = reader.read(wkt); //����geom
		geom->setSRID(indexb); //����S��geom��SRID	
		const Envelope* envelope = geom->getEnvelopeInternal(); //����MBR
		vector<void*> result;
		treea.query(envelope, result); //����S��MBR��R��R�����в�ѯ		
		for (int i = 0; i < result.size(); i++)
		{
			Geometry* geoma = (Geometry*)result[i];
			const Envelope* envelopea = geoma->getEnvelopeInternal(); //����result��MBR 
			Envelope intersect;
			bool flag = envelope->intersection(*envelopea, intersect); //MBR�Ľ���			
			if (flag)
			{
				candidate temp;
				temp.indexa = geoma->getSRID(); //ri
				temp.indexb = indexb; //sj
				temp.centerx = (intersect.getMaxX() + intersect.getMinX()) / 2.0; //pij.x
				temp.centery = (intersect.getMaxY() + intersect.getMinY()) / 2.0; //pij.y
				size_t m = geoma->getNumPoints();
				size_t n = geom->getNumPoints();
				temp.weight = (n + m) * log(n + m); //wij
				C.push_back(temp);
			}
		}
		b.push_back((Geometry*)geom.release());
		indexb++;
	}
	cout << (double)(clock() - astart) / CLOCKS_PER_SEC << endl;
}
void intersect() //���ж����к�ѡ�Խ���intersects��intersection���� 
{
	double sum1, sum2;
	sum1 = sum2 = 0;
	for (int i = 0; i < C.size(); i++)
	{
		astart = clock();
		a[C[i].indexa]->intersects(b[C[i].indexb]);
		sum1 += (double)(clock() - astart) / CLOCKS_PER_SEC;
		astart = clock();
		a[C[i].indexa]->intersection(b[C[i].indexb]);
		sum2 += (double)(clock() - astart) / CLOCKS_PER_SEC;
	}
	cout << sum1 << endl << sum2 << endl;
}
int main()
{
	find_candidates();
	intersect();
	return 0;
}