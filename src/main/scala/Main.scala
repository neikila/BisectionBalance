import cern.colt.matrix.linalg.EigenvalueDecomposition
import cern.colt.matrix.{DoubleFactory2D, DoubleMatrix2D}

/**
  * Created by k.neyman on 10.03.2017.
  */

trait MatrixRepo {
  type Id = Int
  type Mapping = List[Node]
  type Matrix = List[List[Double]]

  type Node = Int
  case class Link(n1: Node, n2: Node)

  def createA: (Mapping, Matrix)
}

trait ExampleMatrix extends MatrixRepo {
  def createA: (Mapping, Matrix) = {
    val matrixARaw =
      (0 :: 0 :: 1 :: 0 :: 1 :: 0 :: Nil) ::
        (0 :: 0 :: 0 :: 1 :: 0 :: 1 :: Nil) ::
        (1 :: 0 :: 0 :: 1 :: 1 :: 0 :: Nil) ::
        (0 :: 1 :: 1 :: 0 :: 0 :: 1 :: Nil) ::
        (1 :: 0 :: 1 :: 0 :: 0 :: 0 :: Nil) ::
        (0 :: 1 :: 0 :: 1 :: 0 :: 0 :: Nil) ::
        Nil
    (List.tabulate(6)(id => id + 1), matrixARaw.map(_.map(_.toDouble)))
  }
}

trait MyVarMatrix extends MatrixRepo {

  private def excluded: List[Node] = 1 :: 2 :: 5 :: 6 ::
    9 :: 10 :: 13 :: 14 :: 15 :: 16 ::
    25 :: 26 :: 27 :: 28 :: 29 :: 30 :: 31 :: 32 :: Nil
  private def extraLinks: List[Link] = Link(8, 21) :: Nil

  private val size = 4

  private def leftPart: (List[Node], List[Link]) = leftPart(excluded)

  private val size2: Node = size * size

  private def leftPart(excluded: List[Node]): (List[Node], List[Link]) = {
    val nodes: List[Node] = (1 to size2).filterNot(excluded contains).toList
    val links = nodes.flatMap { node =>
      val top: Option[Node] = if (node - size > 0) Some(node - size) else None
      val bottom: Option[Node] = if (node + size < size2 + 1) Some(node + size) else None
      val left: Option[Node] = if ((node - 1) % size != 0) Some(node - 1) else None
      val right: Option[Node] = if ((node + 1) % size != 1) Some(node + 1) else None
      (left :: right :: top :: bottom :: Nil).flatten.filterNot(excluded contains).map(Link(node, _))
    }
    (nodes, links)
  }

  private def rightPart: (List[Node], List[Link]) = {
    val (nodes, links) = leftPart(excluded.filter(_ > size2).map(_ - size2))
    (nodes.map(_ + size2), links.map { case Link(n1, n2) => Link(n1 + size2, n2 + size2) })
  }

  private def withExtra: (List[Node], List[Link]) = {
    val (lNodes, lLinks) = leftPart
    val (rNodes, rLinks) = rightPart
    ((lNodes ::: rNodes).sorted, lLinks ::: rLinks :::
        extraLinks ::: extraLinks.map { case Link(n1, n2) => Link(n2, n1) })
  }

  def createA: (Mapping, Matrix) = {
    val (nodes, links) = withExtra
    val groupBy: Map[Node, List[Link]] = links.groupBy(_.n1)

    def createRow(id: Int): List[Double] = {
      val curNode = nodes(id)
      val linksOfNode = groupBy(curNode)
      List.tabulate(nodes.size) { id2 =>
        if (linksOfNode.contains(Link(curNode, nodes(id2)))) 1 else 0
      }
    }

    val matrix: Matrix = List.tabulate(nodes.size)(createRow)
    (nodes, matrix)
  }
}

trait LMatrixCreator extends MatrixRepo {
  var mapping: Mapping = _

  def createLAsArray: Array[Array[Double]] = {
    val (mapping, matrixA) = createA
    this.mapping = mapping
    val matrixBRaw = matrixA.map(_.sum)
    convert(matrixA, matrixBRaw).map(_.toArray).toArray
  }

  def convertRow(row: List[Double], rowNum: Int, diagonal: Double) = {
    row.zipWithIndex.map {
      case (value, i) if i == rowNum => diagonal - value
      case (value, _) => -value
    }
  }

  def convert(m: List[List[Double]], bVector: List[Double]): List[List[Double]] = convert(m, 0, bVector)

  def convert(m: List[List[Double]], rowNum: Int, bVector: List[Double]): List[List[Double]] = {
    m match {
      case h :: t => convertRow(h, rowNum, bVector.head) :: convert(t, rowNum + 1, bVector.tail)
      case _ => Nil
    }
  }
}

object Main extends LMatrixCreator with MyVarMatrix {

  def main(args: Array[String]): Unit = {
    val factory = DoubleFactory2D.dense

    val matrix2D: DoubleMatrix2D = factory.make(createLAsArray)
    println("A")
    var i = 0
    var j = 0
    matrix2D.toArray.foreach { line =>
      j = 0
      line.foreach { el =>
        print(f"${if (i == j) 0 else -el.toInt}%d;")
        j += 1
      }
      i += 1
      println
    }
    println("B")
    i = 0
    j = 0
    matrix2D.toArray.foreach { line =>
      j = 0
      line.foreach { el =>
        print(f"${if (i != j) 0 else el.toInt}%d;")
        j += 1
      }
      i += 1
      println
    }
    println("L")
    matrix2D.toArray.foreach{ line => line.foreach(el => print(f"${el.toInt}%d;")); println }

    val solver = new EigenvalueDecomposition(matrix2D)
    val eigenVectors: DoubleMatrix2D = solver.getV

    mapping.foreach(el => print(s"Q$el;"))
    println
    println("Eigen vals")
    solver.getRealEigenvalues.toArray.foreach(el => print(f"$el%1.2f;"))
    println
    println("Eigen vectors")
    eigenVectors.toArray.foreach{ line => line.foreach(el => print(f"$el%1.2f;")); println }

    val uList: List[(Double, Int)] = eigenVectors.viewColumn(1).toArray.toList.zipWithIndex
    val average: Double = uList.map(_._1).average
    val (uListMedian, uListWithoutMedian) = uList.partition { case (value, id) => value == average }

    val (firstPart, secondPart) = uListWithoutMedian.partition { case (value, id) => value > average }
    println(s"First: ${firstPart.map { case (value, id) => mapping(id)}}")
    println(s"Second: ${secondPart.map { case (value, id) => mapping(id)}}")
    println(s"Middle: ${uListMedian.map { case (value, id) => mapping(id)}}")

    println(s"average = ${average * math.pow(eigenVectors.size(), 0.5)}")
  }

  implicit class ListExtended(list: List[Double]) {
    def average = list.sum / list.size
  }
}
